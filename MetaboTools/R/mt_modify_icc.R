#' mt_modify_icc
#'
#' Calculate ICC - intraclass correlation
#'
#' @param D \code{SummarizedExperiment} input
#' @param qc_samples Logical expression. Can use fields from \code{colData()}
#' @param col_lab Label for the new icc values to be stored in in rowData
#' @param grp_col it is the name of the column in colData containing groups of qc_samples
#' @param icc_lmer Optional flag to set lmer flag in ICC function default is FALSE
#' 
#' @return rowData(D)[[col_lab]] contains column with ICC score for each metabolite
#' @examples
#' \dontrun{... %>% mt_modify_icc(qc_samples=="PQC", col_lab = "PQC_cv", grp_col='dup_grp') %>% ...}
#'
#' @author RB
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#' @importFrom psych ICC
#' @importFrom utils stack
#' @export
mt_modify_icc <- function(D, qc_samples, col_lab, grp_col, icc_lmer=F){
  
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(qc_samples)) stop("qc_samples can't be empty!")
  if(missing(grp_col)) stop("grp_col can't be empty!")
  if(missing(col_lab)) stop("col_lab can't be empty!")
  
  ## APPLY FILTER TO ROW DATA
  qc_samples_q <- dplyr::enquo(qc_samples)
  cd <- colData(D) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("colnames") %>%
    dplyr::filter(!!qc_samples_q)
  
  ## SUBSET SAMPLES
  D1 <- D[, as.numeric(as.matrix(cd$colnames))]
  
  # get the assay of subset with grp column
  qc_data <- bind_cols(D1 %>% assay() %>% t() %>% data.frame(),
                       D1 %>% colData() %>% data.frame() %>% dplyr::select(!!grp_col) %>% data.frame())
  
  # names of the metabolites
  mets <- D1 %>% assay() %>% t() %>% data.frame() %>% colnames()
  
  # ICC computation per metabolite
  icc_df <- lapply(mets, function(x){
    icc_met <- qc_data %>% dplyr::select_at(vars(starts_with(!!grp_col), x))
    f <- as.formula(glue::glue(x, ' ~', sym(grp_col)))
    icc_met <- utils::unstack(icc_met, f)
    icc_val <- tryCatch(psych::ICC(na.omit(icc_met), lmer=F)$results[3, 2], silent=T,error = function(err){NA})
    return(icc_val)
  }) %>% do.call(rbind,.) %>% data.frame(icc=.)
  
  # output in rowData
  rowData(D)[[col_lab]] <- icc_df$icc
  
  ## add status information
  funargs <- MetabolTools:::mti_funargs()
  metadata(D)$results %<>%
    MetabolTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Added QC ICC to rowData!")
    )
  ## return
  D
}
