#' Scale each sample relative to the mean of a given group of samples.
#' 
#' Usually used to make all concentrations relative to a control group, for example.
#' If data are logged, will use a 'minus' operation, if data are non-logged, will divide.
#' Whether data are logged is determine by whether they have been logged inside the MT pipeline.
#'
#' @param D \code{SummarizedExperiment} input
#' @param refsamples expression filtering reference samples from colData
#' @param islogged data logged?
#'
#' @return assay: relatively scaled data
#'
#' @examples
#' # normalize to control group, data not logged
#' ... %>% mt_pre_trans_relative(refsamples = GROUP=="ctrl", islogged=F) %>% ...    
#' 
#' @author JK
#' 
#' @export
mt_pre_trans_relative <- function(
  D,        # SummarizedExperiment input
  refsamples,
  islogged
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  X = t(assay(D))
  
  # check if data are logged
   if (!islogged && (any(unlist(X)[!is.na(unlist(X))]<0))) stop("Data is not logged but contains negative values.")
  
  # find samples to normalize to
  sample_filter_q <- enquo(refsamples)
  cd <- colData(D) %>%
    as.data.frame() %>%
    rownames_to_column("colnames") %>%
    filter(!!sample_filter_q)
  # define samples to be used
  useref <- (colnames(D) %in% cd$colnames) 
  
  # pick operation
  if (!islogged) op <- `/`
  else op = `-`
  # calculate average vector across reference samples
  avg <- X[useref,] %>% apply(2, mean, na.rm=T)
  X %<>% apply(1, function(s)op(s,avg)) # comes out transposed the right way
  
  # save assay
  assay(D) = X
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue('scaled samples relative to {sum(useref)} reference samples: {enquo(refsamples) %>% as.character()}')
    )
  
  # return
  D
  
}
