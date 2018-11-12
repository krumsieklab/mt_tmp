# MetaboTools
#
# Quotient implementation 
# according to Dieterle et al., 2006; https://www.ncbi.nlm.nih.gov/pubmed/16808434
#
# last update: 2018-10-11
# authors: JK
#

require(glue)

mt_pre_norm_quot = function(
  D,                 # SummarizedExperiment input
  vars=1:dim(D)[1],  # vars: index vector of variables of be used, default: all
  NAerror=F,         # NAerrors: throw error for NA's or just ignore?
  refsamples=NULL    # expression filtering reference samples from colData
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
  
  # apply reference sample filter?
  if (!missing(refsamples)) {
    sample_filter_q <- enquo(refsamples)
    cd <- colData(D) %>%
      as.data.frame() %>%
      rownames_to_column("colnames") %>%
      filter(!!sample_filter_q)
    # define samples to be used
    useref <- (colnames(D) %in% cd$colnames) 
  } else {
    useref = rep(T, ncol(D))
  }
  
  # median reference sample
  ref = apply(X[useref,vars],2,function(x)median(x,na.rm=T))
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
      logtxt = glue('quotient normalization based on {sum(useref)} reference samples'),
      output = list(dilution=d)
    )
  
  # return
  assay(D) = t(Y)
  D
  
}
