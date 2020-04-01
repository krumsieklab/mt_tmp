library(glue)

#' Quotient normalization
#'
#' Implementation according to Dieterle et al., 2006\cr
#' https://www.ncbi.nlm.nih.gov/pubmed/16808434
#'
#' @param D \code{SummarizedExperiment} input
#' @param vars index vector of variables of be used, default: all
#' @param NAerror T/F, throw error for NA's or just ignore?
#' @param refsamples expression filtering reference samples from colData
#' 
#' @return assay: quotient-normalized version
#' @return $result: dilution factors
#'
#' @examples
#' #' # in the context of a SE pipeline
#' ... %>% mt_pre_norm_quot() %>% ...    # standard call
#' ... %>% mt_pre_norm_quot(refsamples = GROUP=="ctrl") %>% ...    # use reference samples where 'GROUP' field in colData is 'ctrl'
#' #' ... %>% mt_pre_norm_quot(metMax = 0.2) %>% ...    # use only metabolites with <= 20% missing values to compute the reference used for normalization
#'
#' @author JK
#' 
mt_pre_norm_quot = function(
  D,                 # 
  vars=1:dim(D)[1],  # vars:
  NAerror=F,         # NAerrors: 
  refsamples=NULL,   # 
  metMax=1           # maximum rate of missingness in order to include metabolite in reference 
) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(!(metMax<0 || metMax>1))
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
    # if no samples are left, throw error
    if (sum(useref) == 0) stop(sprintf("No samples match filter for quotient normalization: %s", quo_name(sample_filter_q)))
  } else {
    useref = rep(T, ncol(D))
  }

  # compute metabolite missingness rate
  metmiss <- sapply(1:dim(X)[2], function(k){
    sum(is.na(X[,k]))/dim(X)[1]
  })
  # filter metabolites with missingness rate greater than metMax
  vars <- vars[metmiss <= metMax]
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
      logtxt = glue('quotient normalization based on {sum(useref)} reference samples and {length(vars)} variables: {enquo(refsamples) %>% as.character()}'),
      output = list(dilution=d)
    )
  
  # return
  assay(D) = t(Y)
  D
  
}
