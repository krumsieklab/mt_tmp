#' Identifies sample outliers
#'
#' Uses either univariate or multivariate (leverage, mahalanobis) approaches, won't do both in one call. 
#' The function supports three modes of operation:
#' 1. univariate: a metabolite value is defined as an outlier if it is more than _thresh_ standard deviations from the mean. A sample is defined an outlier if more than _perc_ of its metabolites are univariate outliers.
#' 2. leverage: multivariate approach that uses leverage to define outliers. A sample is defined an outlier if its leverage is larger than _thresh_ times m/n (where m is the number of metabolites and n is the number of samples in the dataset).
#' 3. mahalanobins: multivariate approach that uses the Mahalanobis distance to define outliers. A sample is defined an outlier if its Mahalanobis distance is in the _pval_ quantile of the chi-square distribution.
#'
#' @param D \code{SummarizedExperiment} input
#' @param method Can be either "univariate", "leverage" or "mahalanobis"
#' @param reduce_dim boolean, if TRUE performs PCA-based dimensionality reduction before outlier detection. Can be used to apply multivariate outlier detection methods to low-rank datasets.
#' @param thresh numeric. For the univariate method, number of standard deviations to use as threshold to define an outlier. For the leverage method, number of m/n units to use as threshold to define an outlier. Default value 4.
#' @param perc numeric. For the univariate method, ratio of metabolites that need to be outliers in order to consider the whole sample an outlier. Default value 0.5.
#' @param pval numeric between 0 and 1. For the mahalanobis method, quantile of the chi-squared distribution to use as threshold to define a sample an outlier. Default value 0.01.
#'
#' @return colData: New columns including a binary vector and a numeric score vector.
#' @return $results$output: Returns the specific parameters used to determine outliers for the method selected.
#'
#' @examples
#' \dontrun{# first identify samples that have more than 50% univariate outliers, then identify multivariate outliers with a leverage >4m/n
#' ... %>%
#'   mt_pre_outlier_detection(method="univariate", thresh=4, perc=0.5) %>%
#'   mt_pre_outlier_detection(method="leverage", thresh=4) %>%
#'   mt_pre_outlier_detection(method="mahalanobis", pval=0.01) %>%
#' ...}
#'
#' @author EB, JK
#'
#' @export
mt_pre_outlier_detection <- function(D, method="univariate", reduce_dim=F, ...) {

  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  if(!(method %in% c("univariate","leverage", "mahalanobis")))
    stop("method can only be either univariate, leverage, or mahalanobis")

  X <- t(assay(D))
  X <- scale(X)

  if(any(is.na(X)) & method %in% c("leverage","mahalanobis"))
    stop("Missing values found in the data matrix. Multivariate outlier detection approaches cannot handle missing values.")

  # perform dimension reduction?
  if (reduce_dim) {
    ## Calculate the number of "independent features"
    ## As in Li and Ji, Heredity, 2005
    cordat <- stats::cor(X)
    eigenvals <- eigen(cordat)$values
    Meff <- sum( as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)) )
    ## reduce
    pca <- stats::prcomp(X)
    X <- pca[]
  }

  args <- list(...)

  if(method=="univariate") {
    if(is.null(args$thresh)) args$thresh <- 4
    if(is.null(args$perc)) args$perc <- 0.5
    # compute univariate outliers
    H <- matrix(0, dim(X)[1], dim(X)[2])
    H[X>=args$thresh] <- 1
    # compute percentage of univariate outliers per sample
    score <- rowSums(H)/dim(H)[2]
    # define outliers
    out <- rep(FALSE,length(score))
    out[score > args$perc] <- TRUE

  } else if (method == "leverage") {

    if(!limma::is.fullrank(X))
      stop("The data matrix is not full-rank, leverage cannot be computed")

    if(is.null(args$thresh)) args$thresh <- 4
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
    out <- rep(FALSE,length(score))
    out[score > args$thresh*sum(score)/dim(X)[1]] <- TRUE

  } else if (method == "mahalanobis") {

    if(!limma::is.fullrank(X))
      stop("The data matrix is not full-rank, Mahalanobis cannot be computed")

    if(is.null(args$pval)) args$pval <- 0.01

    # Calculate covariance matrix
    cov_mat <- stats::cov(X)
    # Get mahalanobis distance
    score <- stats::mahalanobis(X, colMeans(X), cov_mat)
    # Define outliers based on threshold
    out <- rep(FALSE,length(score))
    out[score > stats::qchisq(1 - args$pval/nrow(X), df=ncol(X))] <- TRUE

  } else {
    stop("BUG")
  }
  # adding to colData
  colData(D)$outlier <- out
  colData(D)$score <- score
  colnames(colData(D))[colnames(colData(D))=="outlier"] <- paste0("outlier_", method, sep="")
  colnames(colData(D))[colnames(colData(D))=="score"] <- paste0("score_", method, sep="")

  l <- list(threshold=args$thresh)
  if(method=="univariate") l$perc <- args$perc
  names(l) <- paste0(names(l), "_", method, sep="")

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("flagged %d %s outliers", sum(out, na.rm = TRUE), method),
      output = l
    )

  # return
  D

}

# deleted code for low-rank cov approx
# is_invertible <- class(try(solve(cov_mat),silent=T))=="matrix"
# if (!is_invertible) {
#   # matrix has low rank, use GeneNet's shrinkage version instead
#   # https://rdrr.io/cran/Hotelling/man/hotelling.stat.html
#   cov_mat <- cov.shrink(X, verbose=F) %>% apply(2, as.vector)
# }
#
#
#
#

#
# ##### SHRINKAGE COV and PINV
#
# # Calculate covariance matrix and inverse
# C <- corpcor::cov.shrink(X, verbose=F) %>% apply(2, as.vector)
# # Get mahalanobis distance
# score <- stats::mahalanobis(X, colMeans(X), C)
#
#
# # # Xc <- scale(X, scale=F, center=T)
# # center = colMeans(X)
# # Xc <- sweep(X, 2L, center)
# # score2 <- rowSums(Xc %*% Ci * Xc)
# # i=5; t(Xc[i,]) %*% Ci %*% Xc[i,]
# # Define outliers based on threshold
# out <- rep(FALSE,length(score))
# out[score > stats::qchisq(1 - args$pval/nrow(X), df=ncol(X))] <- TRUE
