#' Replace SummarizedExperiment with another one.
#' 
#' This function replaces the assay(), colData() and rowData() of the current pipeline object with those of another one. 
#' The metadata of the original pipeline object is retained.
#' 
#' The function provides a possibility for "pipeline restarts" after operations that entirely change the data matrix, such
#' as mt_modify_ratios and mt_modify_aggPW.
#'
#' @param D SummarizedExperiment object
#' @param Drepl replacement SummarizedExperiment object where assay(), colData() and rowData() will be copied from.
#'
#' @examples
#' \dontrun{# Transform dataset to ratios
#' ... %>%  mt_modify_replaceSE(Drepl = D1) %>% ... # replace with an earlier SE
#' }
#'
#' @return SummarizedExperiment containing assay(), colData() and rowData() of Drepl, but metadata of the original pipeline object
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
#' 
mt_modify_replaceSE <- function(D, Drepl) {
  
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot("SummarizedExperiment" %in% class(Drepl))
  
  ## CREATE NEW OBJECT
  Dnew <- SummarizedExperiment(
    assay    = assay(Drepl),
    rowData  = rowData(Drepl),
    colData  = colData(Drepl),
    metadata = metadata(D))
  
  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = "replaced SE object",
      output = NULL
    )
  Dnew
  
  
}
