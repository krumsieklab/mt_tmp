# MetaboTools
#
# Exponentiate, base 2 by default.
#
# last update: 2018-10-16
# authors: JK
#

mt_pre_trans_exp <- function(
  D,      # SummarizedExperiment input
  base=2  # base of exp
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(base%%1==0) # integer number
  
  # exp
  assay(D) = base^(assay(D))
  
  # add status information
  call = match.call()
  metadata(D)$preprocess %<>% 
    add_to_list(list(txt=sprintf('exp, base %d', base), call=call))
  
  # return
  D
  
}
