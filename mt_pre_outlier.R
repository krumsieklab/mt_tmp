require(limma)
source(codes.makepath("MT/mt_internal_helpers.R"))

#' Identifies sample outliers.
#' 
#' Uses either univariate or multivariate (leverage) approaches, won't do both in one call.
#'
#' @param D \code{SummarizedExperiment} input
#' @param method Can be either "univariate" or "leverage" (for now)
#' @param thr Number of standard deviations or m/n units to use as threshold to define the outlier
#' @param perc For the univariate method, percentage of metabolites that need to be outliers in order to consider the whole sample an outlier
#'
#' @return SE rows or columns will be filtered.
#' @return $output: logical vector of outlier samples as well as a numeric vector with the sample scores
#'
#' @examples
#' # first identify samples that have more than 50% univariate outliers, then identify multivariate outliers with a leverage >4m/n
#' ... %>%
#'   mt_pre_outlier(method="univariate", thr=4, perc=0.5) %>%
#'   mt_pre_outlier(method="leverage", thr=4) %>%
#' ...
#' 
#' @author EB
#' 

mt_pre_outlier <- function(
  D,            # SummarizedExperiment input
  method="leverage",     # Method for outlier detection, can only be either "univariate" or "leverage"
  thr=4,        # Number of standard deviations (univariate) or of m/n units (leverage) to use as threshold for the definition of outlier
  perc = 0.5    # For the univariate method, percentage of metabolites that need to be outliers in order to consider the whole sample an outlier
) {
  
  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  if(!(method %in% c("univariate","leverage")))
    stop("method can only be either univariate or leverage")
  
  X <- t(assay(D))
  X <- scale(X)
  
  if(any(is.na(X)))
    stop("Missing values found in the data matrix")
  if(method=="leverage" & !is.fullrank(X))
    stop("The data matrix is not full-rank, leverage cannot be computed")

  if(method=="univariate") {
    # compute univariate outliers
    H <- matrix(0, dim(X)[1], dim(X)[2])
    H[X>=thr] <- 1
    # compute percentage of univariate outliers per sample
    score <- rowSums(H)/dim(H)[2]
    # define outliers
    out <- rep(0,length(score))
    out[score > perc] <- 1
  } else {
    # # compute hat matrix through Singular Value Decomposition
    # SVD <- svd(X)
    # H <- tcrossprod(SVD$u)
    
    # compute hat matrix with the Pivoted Cholesky factorization
    L <- t(suppressWarnings(chol(crossprod(X), pivot = TRUE)))
    r <- attr(L, "rank")
    piv <- attr(L, "pivot")
    Qt <- forwardsolve(L, t(X[, piv]), r)
    H <- crossprod(Qt)
    
    # extract leverage values
    score <- diag(H)
    # define outliers
    out <- rep(0,length(score))
    out[score > thr*sum(score)/dim(X)[1]] <- 1 
  }
  
  # adding to colData
  colData(D)$outlier <- out
  colData(D)$score <- score
  colnames(colData(D))[colnames(colData(D))=="outlier"] <- paste0("outlier_", method, sep="")
  colnames(colData(D))[colnames(colData(D))=="score"] <- paste0("score_", method, sep="")
  
  l <- list(threshold=thr)
  if(method=="univariate") l$perc <- perc
  names(l) <- paste0(names(l), "_", method, sep="")
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("flagged %d %s outliers", sum(out), method),
      output = l
    )
  
  # return
  D
  
}