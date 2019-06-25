require(limma)
source(codes.makepath("MT/mt_internal_helpers.R"))

#' Identifies sample outliers.
#' 
#' Uses either univariate or multivariate (leverage) approaches, won't do both in one call.
#'
#' @param D \code{SummarizedExperiment} input
#' @param threshold Number of standard deviations or m/n units to use as threshold to define the outlier, default value set to 4
#'
#' @return SE with NA values where outliers used to be
#'
#' @examples
#' # first identify samples that have more than 50% univariate outliers, then identify multivariate outliers with a leverage >4m/n
#' ... %>%
#'   mt_pre_singleout(threshold=3) %>%
#' ...
#' 
#' @author Annalise Schweickart
#' 

mt_pre_outliercorrection <- function(
  D,            # SummarizedExperiment input
  threshold=4,  # threshold for outlier detection (default is 4 standard deviations)
  ...   
) {
  
  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  
  X <- t(assay(D))
  X <- scale(X)
  
  if(any(is.na(X)))
    stop("Missing values found in the data matrix")

  # compute univariate outliers
  H <- matrix(F, dim(X)[1], dim(X)[2])
  H[abs(X)>=threshold] <- T
  # change outliers to NA
  assay(D)[t(H)] <- NA
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("flagged %d outliers", sum(H))
    )
  
  # return
  D
  
}