#### Reorder a sample annotation
# will convert it to a factor first

library(gdata) 

mt_modify_sampleanno_reorder <- function(
  D,          # SummarizedExperiment input
  field,      # field to access
  new.order   # new order, either numeric vector with new order, or has to contain all fields
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # get variable
  if (!(field %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", field))
  p = colData(D)[[field]]
  
  # ensure factor
  p <- as.factor(p)
  
  # different actions for numeric or factor orders
  if (is.numeric(new.order)) {
    if (length(new.order) != length(levels(p)) || !(all.equal(sort(new.order), 1:length(levels(p)))==T)) stop(sprintf("For numeric new.order, the vector has to contain all numbers from 1 to number of levels. Expected: 1:%d",length(levels(p))))
    p <- gdata::reorder.factor(p, new.order=new.order)
  } else {
    if (length(new.order) != length(levels(p)) || !(all.equal(sort(new.order), sort(levels(p)))==T)) stop(sprintf("For string new.order, the vector has to contain all levels exactly once. Levels: %s", paste0(levels(p),collapse=',')))
    p <- gdata::reorder.factor(p, new.order=new.order)
  }
  
  # write back
  colData(D)[[field]] <- p
  
  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Reordered column '%s' as '%s'", field, paste0(new.order,collapse=','))
    )
  ## return
  D
  
}






