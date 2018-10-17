# MetaboTools
#
# Simple logging of data, base 2 by default.
#
# last update: 2018-10-12
# authors: JK
#

mt_pre_trans_log <- function(
  D,      # SummarizedExperiment input
  base=2  # base of logging
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(base%%1==0) # integer number
  
  # log
  assay(D) = log(assay(D))
  
  # add status information
  call = match.call()
  metadata(D)$preprocess %<>% 
    add_to_list(list(
      type='log',
      label=sprintf("log%d",base),
      log=sprintf("log%d",base), 
      call=call
      ))
  
  # return
  D
  
}
