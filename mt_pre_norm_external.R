library(glue)

#' Normalize by an external sample annotation
#' 
#' Normalize by a field in colData, e.g. sample weight, DNA content, protein content (BRADFORD) etc.
#'
#'
#' @param D \code{SummarizedExperiment} input
#' @param field numeric field in \code{colData} (samples) to normalize by
#' 
#' @return assay: normalized version
#'
#' @examples
#' #' # in the context of a SE pipeline
#' ... %>% mt_pre_norm_quot(field='DNA') %>% ...    # normalize by values in field DNA
#'
#' @author JK
#' @export
#' 
mt_pre_norm_external = function(
  D,                 # 
  field
) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))
  if (any(unlist(X)[!is.na(unlist(X))]<0)) stop("Matrix contains negative values. Did you input logged data?")
  
  # get variable to normalize by
  if (!(field %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", field))
  vc = colData(D)[[field]]
  if (!is.numeric(vc)) stop(sprintf("'%s' is not numeric.", field))
  
  # run normalization
  Y = t(sapply(1:dim(X)[1], function(i)unlist(X[i,]/vc[i])))
  rownames(Y) = rownames(X)
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue("normalized by '{field}'")
    )
  
  # return
  assay(D) = t(Y)
  D
  
}
