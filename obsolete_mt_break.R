# MetaboTools
#
# This will add a non-action entry with a unique code to the pipeline protocol.
# Should be used to merge results from multiple pipeline arms later.
#
# last update: 2018-11-12
# authors: JK
#

library(uuid)
library(glue)

mt_break <- function(
  D
) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # add status information & plot
  id <- UUIDgenerate()
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue("pipeline break, {id}"),
      output = id
    )
  
  # return
  D
  
}


