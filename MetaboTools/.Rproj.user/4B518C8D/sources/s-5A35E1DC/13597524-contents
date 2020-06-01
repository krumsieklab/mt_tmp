#' Identifies sample outliers.
#'
#' Uses either univariate or multivariate (leverage) approaches, won't do both in one call.
#'
#' @param D \code{SummarizedExperiment} input
#' @param method Can be either "univariate" or "leverage" (for now)
#' @param reduce_dim Perform PCA-based dimension reduction before outlier detection? Needed for multivariate methods in low-rank datasets.
#' @param thresh Number of standard deviations or m/n units to use as threshold to define the outlier, default value set to 4
#' @param perc For the univariate method, percentage of metabolites that need to be outliers in order to consider the whole sample an outlier, default value set to 0.5
#' @param pval For the mahalanobis distance method, p-val of chi-squared test to threshold at, default = 0.01
#'
#' @return SE with additional colData columns including a binary vector and a numeric score vector
#' @return $output: returns the specific parameters used to determine outliers for the specific method selected
#'
#' @examples
#' \dontrun{# first identify samples that have more than 50% univariate outliers, then identify multivariate outliers with a leverage >4m/n
#' ... %>%
#'   mt_pre_outlier(method="univariate", thresh=4, perc=0.5) %>%
#'   mt_pre_outlier(method="leverage", thresh=4) %>%
#'   mt_pre_outlier(method="mahalanobis", pval=0.01) %>%
#' ...}
#'
#' @author EB, JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_pre_outlier <- function(
  D,            # SummarizedExperiment input
  method="univariate",     # Method for outlier detection, can only be either "univariate" or "leverage"
  reduce_dim=F,
  ...
) {

  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  if(!(method %in% c("univariate","leverage", "mahalanobis")))
    stop("method can only be either univariate, leverage, or mahalanobis")

  X <- t(assay(D))
  X <- scale(X)

  if(any(is.na(X)))
    stop("Missing values found in the data matrix")

  # perform dimension reduction?
  if (reduce_dim) {
    ## Calculate the number of "independent features"
    ## As in Li and Ji, Heredity, 2005
    cordat <- cor(X)
    eigenvals <- eigen(cordat)$values
    Meff <- sum( as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)) )
    ## reduce
    pca <- prcomp(X)
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
    cov_mat <- cov(X)
    # Get mahalanobis distance
    score <- mahalanobis(X, colMeans(X), cov_mat)
    # Define outliers based on threshold
    out <- rep(FALSE,length(score))
    out[score > qchisq(1 - args$pval/nrow(X), df=ncol(X))] <- TRUE

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
# score <- mahalanobis(X, colMeans(X), C)
#
#
# # # Xc <- scale(X, scale=F, center=T)
# # center = colMeans(X)
# # Xc <- sweep(X, 2L, center)
# # score2 <- rowSums(Xc %*% Ci * Xc)
# # i=5; t(Xc[i,]) %*% Ci %*% Xc[i,]
# # Define outliers based on threshold
# out <- rep(FALSE,length(score))
# out[score > qchisq(1 - args$pval/nrow(X), df=ncol(X))] <- TRUE
