#### store info to add heading during report generation

mt_reporting_heading <- function(
  D,
  strtitle,
  lvl=1
) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # add status information & heading info
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue("reporting heading, level {lvl}: {strtitle}"),
      output = list(lvl=lvl,title=strtitle)
    )
  
  # return
  D
  
}


