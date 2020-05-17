#' Log, base 2 by default.
#' 
#' Transform the entire dataset log_based(x), i.e. log(x)/log(base)
#'
#' @param D \code{SummarizedExperiment} input
#' @param base operation: log(x)/log(base) for every data point
#'
#' @return assay: logged data
#'
#' @examples
#' # in the context of a SE pipeline
#' ... %>% mt_pre_trans_log() %>% ...    # standard call, base 2
#' ... %>% mt_pre_trans_log(base=10) %>% ...    # base 10
#' 
#' @author JK
#' 
mt_pre_trans_log <- function(
  D,      # SummarizedExperiment input
  base=2  # base of logging
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(base%%1==0) # integer number
  
  # log
  assay(D) = log(assay(D), base = base)
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("log%d", base)
    )
  
  # return
  D
  
}
