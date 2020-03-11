require(limma)
source(codes.makepath("MT/mt_internal_helpers.R"))

#' Identifies single outliers in samples.
#' 
#' Uses univariate approach
#'
#' @param D \code{SummarizedExperiment} input
#' @param threshold Number of standard deviations or m/n units to use as threshold to define the outlier, default value set to 4
#' @param sample_num_correction Whether number of outliers should depend on number of samples or hard sd cutoff
#'
#' @return SE with NA values where outliers used to be
#'
#' @examples
#' ... %>%
#'   mt_pre_outliercorrection(threshold=3) %>%
#' ...
#' 
#' @author Annalise Schweickart
#' 

mt_pre_outliercorrection <- function(
  D,            # SummarizedExperiment input
  threshold=NA,  # threshold for outlier detection (default is 4 standard deviations)
  sample_num_correction = T, # Should the number of outliers be corrected by the sample size
  ...   
) {
  
  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  
  X <- t(assay(D))
  X <- scale(X)
  # 
  # if(any(is.na(X)))
  #   stop("Missing values found in the data matrix")
  if(is.na(threshold) & sample_num_correction == F){
    stop("Threshold must be provided if not corrected by sample numbers")
  }
  if(is.na(threshold)){
    numsamp=nrow(X)
    alpha=0.05/numsamp
    threshold = qnorm( 1 - (alpha/2) )
  }
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