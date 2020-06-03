#' mt_modify cv
#'
#' Calculate CV score
#'
#' @param D \code{SummarizedExperiment} input
#' @param qc_samples Logical expression. Can use fields from \code{colData()}
#' @param col_lab Label for the new cv values to be stored in in colData
#'
#' @return rowData(D)[[col_lab]] contains column with CV score for each metabolite
#' @examples
#' \dontrun{... %>% mt_modify_cv(qc_samples=="PQC", col_lab = "PQC_cv") %>% ...}
#'
#' @author Annalise Schweickart
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_modify_cv <- function(D, qc_samples, col_lab){

  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(qc_samples))
    stop("qc_samples can't be empty")

  ## APPLY FILTER TO ROW DATA
  qc_samples_q <- dplyr::enquo(qc_samples)
  cd <- colData(D) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("colnames") %>%
    dplyr::filter(!!qc_samples_q)

    ## SUBSET SAMPLES
  D_cv <- D[, cd$colnames]

  ## Calculate metabolite cv scores
  rowData(D)[[col_lab]] <- apply(assay(D_cv), 1, function(x) stats::sd(x, na.rm=T)/mean(x, na.rm=T))

  ## add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Add QC cv to colData ")
      )
  ## return
  D
}
