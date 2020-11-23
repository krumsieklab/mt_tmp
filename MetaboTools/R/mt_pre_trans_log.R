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
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_trans_log() %>% ...    # standard call, base 2
#' ... %>% mt_pre_trans_log(base=10) %>% ...    # base 10
#' }
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_pre_trans_log <- function(
  D,      # SummarizedExperiment input
  base=2  # base of logging
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if(base%%1!=0) warning(sprintf("Base not integer. Please double check. Base %f.", base))

  # log
  assay(D) = log(assay(D), base = base)

  # add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("log%f", base)
    )

  # return
  D

}
