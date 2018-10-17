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
  call = match.call()
  metadata(D)$preprocess %<>% 
    add_to_list(list(txt=sprintf('scaled, center=%d, scale=%d', center, scale), call=call))
  
  # return
  D
  
}
