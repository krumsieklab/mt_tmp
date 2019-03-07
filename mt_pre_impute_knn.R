# MetaboTools
#
# Impute using kNN method.
#
# Default settings are from the winning approach of our paper.
# https://link.springer.com/article/10.1007%2Fs11306-018-1420-2
#
# last update: 2018-10-12
# authors: JK
#

source(codes.makepath("MT/mt_internal_helpers.R"))
source(codes.makepath("MT/mt_internal_imputeKNN.R"))

mt_pre_impute_knn <- function(D, methods="knn.obs.euc.sel", K=10, verbose=F) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(K%%1==0) # integer number
  
  # impute
  assay(D) = imputeKNN(assay(D), methods=methods, K=K, verbose=verbose)
  
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