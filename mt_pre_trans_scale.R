#' Scale data, mean 0 / sd 1 by default.
#'
#' @param D \code{SummarizedExperiment} input
#' @param center T/F, mean-center data? default: T
#' @param scale T/F, scale data to sd 1? default: T
#'
#' @return assay: scaled data
#'
#' @examples
#' # in the context of a SE pipeline
#' ... %>% mt_pre_trans_scale() %>% ...    # standard call, center and scale
#' ... %>% mt_pre_trans_scale(scale=F) %>% ...    # only mean centering
#' 
#' @author JK
#' 
#' @export
mt_pre_trans_scale <- function(
  D,        # SummarizedExperiment input
  center=T, # mean 0?
  scale=T   # SD 1?
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.logical(center)) 
  stopifnot(is.logical(scale)) 
  
  # scale
  assay(D) = t(scale(t(assay(D)),center=center,scale=scale))
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('scaled, center=%d, scale=%d', center, scale),
      logshort = 'scaled'
    )

  # return
  D
  
}
