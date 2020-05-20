# MT compatible function --------------------------------------------------

#' Impute using kNN method with multi-core functionality.
#'
#' 
#' Default settings are from the winning approach of Kiki's paper paper.
#' https://link.springer.com/article/10.1007%2Fs11306-018-1420-2
#' Specifically "knn.obs.euc.sel"
#' 
#' This script creates a tmp folder where distance matrices and imputed
#' values are stored
#
#' @param D \code{SummarizedExperiment} input
#' @param K Number of nearest neighbors to consider, default is 10 (recommendation: do not touch)
#' @param mc_cores Number of cores to use for imputation, dafualt: 5
#' @param verbose T/F, whether to output intermediate steps, default: F
#'
#' @return assay: imputed data
#'
#' @examples
#' # in the context of a SE pipeline
#' ... %>% mt_pre_impute_knn_multicore() %>% ...    # standard call
#' 
#' @author Parviz Gomari
#' 
#' @export

mt_pre_impute_knn_multicore <- function(D, K=10, mc_cores = 5, verbose=F) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(K%%1==0) # integer number
  
  # impute
  assay(D) =  
    assay(D) %>% 
    t() %>% 
    mt_internal_imputeKNN_multicore(K=K, mc_cores = mc_cores, verbose=verbose) %>% 
    t()

  
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
