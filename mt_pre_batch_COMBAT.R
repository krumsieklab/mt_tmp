# MetaboTools
#
# COMBAT batch correction.
#
# last update: 2018-10-29
# authors: JK
#

require(sva)


mt_pre_batch_COMBAT = function(
  D,       # SummarizedExperiment input
  batches  # sample annotation column that contains batch info
) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = assay(D) # combat uses data probe s sample
  
  # get variable
  if (!(batches %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", batches))
  b = colData(D)[[batches]]
  # ensure it's a factor or character vector
  if (!(is.character(b) || is.factor(b))) stop(sprintf("'%s' has to be character, factor, or numeric.", batches))
  b = as.factor(b)
  
  # run ComBat
  # invisible (capture.output(     # seems to be the only way to silence ComBat() output to stdout and stderr
  #   invisible( capture.output(
      Y <- ComBat(
        dat=X,
        batch=b,
        mod=NULL,
        par.prior=TRUE,
        prior.plots=FALSE)
  #     ), type = c("output")
  #   ), type = c("message")
  # )
  # )
  # )
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("COMBAT batch correction, variable: '%s'",batches)
    )
  
  # return
  assay(D) = Y
  D
  
  
}
