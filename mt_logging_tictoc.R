# MetaboTools
#
# tic/toc functionality to time certain branches of the pipeline.
#
# last update: 2019-01-06
# authors: JK
#

library(tictoc)

mt_logging_tic <- function(D) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # tic
  tic()
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("tic")
    )
  
  # return
  D
  
}

mt_logging_toc <- function(D) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # toc 
  t <- toc(quiet=T)
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("toc, elapsed: %.2fs", t$toc-t$tic)
    )
  
  # return
  D
  
}