# NONFUNCTIONAL
# MetaboTools
#
# Median batch correction, as Metabolon does it for runday correction.
#
# last update: 2018-10-16
# authors: JK
#

mt_pre_batch_median = function(
  D,  # SummarizedExperiment input
  batchvar
) {
  
  stop("NONFUNCTIONAL right now")
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))

  # no negative values allowed
  if (any(unlist(X)[!is.na(unlist(X))]<0)) stop("Matrix contains negative values. Did you input logged data?")
  # check if data actually have been logged by preprocessing
  
  
  # check if batches given, and data unlogged
  if (length(levels(as.factor(batches)))==1) warning("no batches given, performing overall median scaling")
  if (logged) stop("cannot perform median-scaling after logging")
  # median per metabolite
  txt <- printstore(txt,verbose,"Median-scaling per batch")
  X = curdata
  for (i in 1:length(levels(batches))) {
    batch = levels(batches)[i]
    # median normalize
    Xnorm = apply(X[batches==batch,],2,function(c){c/median(c,na.rm=T)})
    X[batches==batch,] = Xnorm
  }
  curdata = X
  
  
}
