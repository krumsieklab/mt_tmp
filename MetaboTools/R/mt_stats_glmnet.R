#' Title / Header of Function
#'
#' Description of function
#'
#' @param D \code{SummarizedExperiment} input
#' @param x
#' @param y
#' @param outer_k
#' @param inner_k
#' @param alpha
#' @param lambda
#'
#' @return what the function returns
#'
#' @examples
#' # example of how to run function
#' \dontrun{}
#'
#' @author KC
#'
#'
#' @noRd

mt_stats_glmnet <- function(D,
                            x,
                            y,
                            outer_k,
                            inner_k,
                            alpha,
                            lambda,
                            rand_seed){

  for(i in 1:outer_k){

    set.seed(rand_seed)
    iMod <- glmnet::cv.glmnet(x=x, y=x, alpha=alpha, lambda=lambda, nfolds=inner_k, family="binomial")

  }








}
