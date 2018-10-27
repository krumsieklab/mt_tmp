# MetaboTools
#
# Quotient implementation 
# according to Dieterle et al., 2006; https://www.ncbi.nlm.nih.gov/pubmed/16808434
#
# last update: 2018-10-11
# authors: JK
#

mt_pre_norm_quot = function(
  D,                 # SummarizedExperiment input
  vars=1:dim(X)[2],  # vars: index vector of variables of be used, default: all
  NAerror=F,         # NAerrors: throw error for NA's or just ignore?
  refsamples=NA      # indices (or logical vector) of samples to calculate reference sample on (e.g. only on control samples)
) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))
  if (any(unlist(X)[!is.na(unlist(X))]<0)) stop("Matrix contains negative values. Did you input logged data?")
  
  # check if there are any NAs
  if (sum(is.na(X[,vars]))>0) {
    # throw warning or error?
    if (NAerror) {
      stop('Data matrix contains NAs')
    } else {
      mti_logwarning('Data matrix contains NAs')
    }
  }
  
  # if reference samples not given -> use all samples
  if (is.na(refsamples)) refsamples <- 1:nrow(X)
  
  # median reference sample
  ref = apply(X[refsamples,vars],2,function(x)median(x,na.rm=T))
  # get dilution factors
  d = apply(X[,vars],1,  function(s) median(as.numeric(s/ref),na.rm=T))
  # apply to each sample  (for each row=sample, divide values by median dilution factor)
  Y = t(sapply(1:dim(X)[1], function(i)unlist(X[i,]/d[i])))
  
  #Y = t(apply(X,1,  function(s) s /  d) )
  rownames(Y) = rownames(X)
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = 'quotient normalization',
      output = list(dilution=d)
    )
  
  # return
  assay(D) = t(Y)
  D
  
}
