#' Impute using kNN method.
#'
#' Default settings are from the winning approach of our paper.
#' https://link.springer.com/article/10.1007%2Fs11306-018-1420-2
#
#' @param D \code{SummarizedExperiment} input
#' @param methods Specific KNN method to use, default is "knn.obs.euc.sel" (recommendation: do not touch)
#' @param K Number of nearest neighbors to consider, default is 10 (recommendation: do not touch)
#' @param verbose T/F, whether to output intermediate steps, default: F
#'
#' @return assay: imputed data
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_impute_knn() %>% ...    # standard call
#' }
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_pre_impute_knn <- function(D, method="knn.obs.euc.sel", K=10, verbose=F) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(K%%1==0) # integer number

  # impute
  assay(D) = t( mti_imputeKNN( t(assay(D)), methods=method, K=K, verbose=verbose) )

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = 'imputed via KNN',
      logshort = 'impute KNN'
    )

  # return
  D

}
