#### Apply function to metabolite or sample annotation column

library(gdata) 

mt_modify_applytoanno <- function(
  D,          # SummarizedExperiment input
  annotype,   # "samples" or "metabolites"
  field,      # field to access
  fun         # function to be applied
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if(!(any(c("samples","metabolites")%in%annotype)))stop("'annotype' must be 'samples' or 'metabolites'")
  
  # get data frame
  df = if(annotype=="samples"){colData(D)}else{rowData(D)}
  
  # get variable
  if (!(field %in% colnames(df))) stop(sprintf("'%s' not found in %s annotations.", field, ifelse(annotype=="samples","sample","metabolite")))
  p = colData(D)[[field]]
  
  # apply function
  pnew <- sapply(p, fun)
  
  # write back
  if (annotype=="samples") {
    colData(D)[[field]] <- pnew
  } else {
    rowData(D)[[field]] <- pnew
  }
  
  
  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Transformed column '%s' of %s annotations", field, ifelse(annotype=="samples","sample","metabolite"))
    )
  ## return
  D
  
}






