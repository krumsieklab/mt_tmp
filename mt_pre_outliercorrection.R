require(limma)
source(codes.makepath("MT/mt_internal_helpers.R"))

#' Identifies single outliers in samples.
#' 
#' Uses univariate approach
#'
#' @param D \code{SummarizedExperiment} input
#' @param threshold Number of standard deviations or m/n units to use as threshold to define the outlier, if not provided, default is correction by sample size. This parameter is only used if sample_num_correction is False
#' @param sample_num_correction Whether number of outliers should depend on number of samples or hard sd cutoff. If true, threshold is ignored
#' @param alpha The percentage of points we will consider outliers based on the assumption of a normal distribution (alpha/2 on either tail). Only used if sample_num_correction is True
#' 
#' @return SE with NA values where outliers used to be
#'
#' @examples
#' ... %>%
#'   mt_pre_outliercorrection(threshold=3, sample_num_correction=F) %>%
#' ...
#' 
#' @author Annalise Schweickart
#' 

mt_pre_outliercorrection <- function(
  D,            # SummarizedExperiment input
  threshold=NA,  # threshold for outlier detection
  sample_num_correction = T, # Should the number of outliers be corrected by the sample size
  alpha = 0.05, # Percent of sample data should be considered an outlier (assuming normal distribution)
  ...   
) {
  
  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  
  X <- t(assay(D))
  X <- scale(X)

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