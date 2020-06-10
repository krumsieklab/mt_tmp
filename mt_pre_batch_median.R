#' Median batch correction.
#'  
#' Same approach that Metabolon uses for runday correction. Divides metabolite values by the median value per batch.
#'
#' @param D \code{SummarizedExperiment} input
#' @param batches # sample annotation (colData) column name that contains batch assignment
#'
#' @return assay: batch-corrected version
#' 
#' @examples
#' ... %>% mt_pre_batch_median(batches="BATCH") %>% ...
#'
#' @author JK
#' 
#' @export
mt_pre_batch_median = function(
  D,       # SummarizedExperiment input
  batches  # sample annotation column that contains batch info
) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))
  
  # get variable
  if (!(batches %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", batches))
  b = colData(D)[[batches]]
  # ensure it's a factor or character vector
  if (!(is.character(b) || is.factor(b))) stop(sprintf("'%s' has to be character, factor, or numeric.", batches))
  b = as.factor(b)
  
  # no negative values allowed
  if (min(X,na.rm=T)<0) stop("Matrix contains negative values.")
  # check if data actually have been logged by preprocessing
  if (length(mti_res_get_path(D, c("pre","trans","log"))) > 0)
    stop("Median batch correction can only be performed on non-logged data.")
  
  # median per metabolite
  for (i in 1:length(levels(b))) {
    batch = levels(b)[i]
    # median normalize
    Xnorm = apply(X[b==batch,,drop=F],2,function(c){c/median(c,na.rm=T)})
    X[b==batch,] = Xnorm
  }

  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("median-scaling per batch in '%s'",batches)
    )
  
  # return
  assay(D) = t(X)
  D
  
  
}
