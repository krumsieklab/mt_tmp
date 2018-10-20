# MetaboTools
#
# Scale data, mean 0, SD 1 by default.
#
# last update: 2018-10-16
# authors: JK
#

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
