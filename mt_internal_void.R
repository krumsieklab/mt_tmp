# function that does not do anything, but leaves an entry in the log
# takes arbitrary arguments and ignores them

mt_internal_void <- function(
  D, param=NA
) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = "void function"
    )
  
  # return
  D
  
}


