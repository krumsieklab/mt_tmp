#' mt_modify_cv
#'
#' Calculate CV - coefficient of variation
#'
#' @param D \code{SummarizedExperiment} input
#' @param qc_samples Logical expression. Can use fields from \code{colData()}
#' @param col_lab Label for the new cv values to be stored in in rowData
#' @param replicates Optional flag to indicate if the selected samples are replicates
#' @param id_col Optional if replicates is T, it is the name of the column in colData containing sample IDs
#' 
#' @return rowData(D)[[col_lab]] contains column with CV score for each metabolite
#' @examples
#' \dontrun{... %>% mt_modify_cv(qc_samples=="PQC", col_lab = "PQC_cv") %>% ...}
#'
#' @author Annalise Schweickart, RB (updated on 2020-09-12)
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_modify_cv <- function(D, qc_samples, col_lab, replicates=F, id_col=NULL){

  stopifnot("SummarizedExperiment" %in% class(D))
  
  if(missing(qc_samples)) stop("qc_samples can't be empty")

  ## APPLY FILTER TO ROW DATA
  qc_samples_q <- dplyr::enquo(qc_samples)
  cd <- colData(D) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("colnames") %>%
    dplyr::filter(!!qc_samples_q)

    ## SUBSET SAMPLES
  D_cv <- D[, as.numeric(as.matrix(cd$colnames))]
  # calc_cv function
  calc_cv <- function(x){sd(x, na.rm = T)/mean(x, na.rm = TRUE)}
  ## Calculate metabolite cv scores
  if (replicates) {
    if(missing(id_col)) stop(sprintf('id_col must be provided if cv is to be calcualted on replicates!'))
    # get the assay of duplicates with id column
    cv_data <- dplyr::bind_cols(D_cv %>% assay() %>% t() %>% data.frame(),
         D_cv %>% colData() %>% data.frame() %>% dplyr::select(!!id_col) %>% data.frame())
    # mean per ID
    rowData(D)[[col_lab]] <- cv_data %>%
      dplyr::group_by_at(vars(starts_with(!!id_col))) %>% # calculate cv per duplicated ID
      dplyr::summarise_at(vars(-starts_with(!!id_col)), calc_cv) %>% select(-!!id_col) %>% 
      colMeans(na.rm = T) %>% unlist()
    
  } else {
    rowData(D)[[col_lab]] <- apply(assay(D_cv), 1, calc_cv)
  }

  ## add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Added QC cv to rowData")
      )
  ## return
  D
}
