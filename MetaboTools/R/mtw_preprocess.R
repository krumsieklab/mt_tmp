#' mtw_preprocess
#'
#' Preprocess on data, including batch correction, quotient normalization, log transformation and kNN imputations
#'
#' @param D \code{SummarizedExperiment} input
#' @param batch_col Batch column if needed
#' @param batch_fun Batch function to use. Default is \code{mt_pre_batch_median}
#' @param boxplot_options A list of parameters for \code{mt_plots_sampleboxplot}
#' @param quot_options A list of parameters for \code{mt_pre_norm_quot}
#' @param comp_fac List of names in colData column to correlation dilution factors with
#' @param dilution_options A list of parameters for \code{mt_plots_qc_dilutionplot}
#' @param log_base Base for log function, default is 2
#' @param do_impute Logical. Whether to run imputation
#' @param knn_options A list of parameters for \code{mt_pre_impute_knn}
#'
#' @return D that has been preprocessed
#' @examples
#' \dontrun{... %>% mtw_preprocess() %>% ...}
#' \dontrun{... %>% mtw_preprocess(batch_col = "BATCH_MOCK",comp_fac = c('num1','Group'), batch_col = "BATCH_MOCK",comp_fac = c("num1","Group"),quot_options = list(ref_samples = quote(Group == 'Vehicle'))) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mtw_preprocess <-
  function(D,
           batch_col,
           batch_fun = mt_pre_batch_median,
           batch_ref_samples = NULL,
           boxplot_options = list(),
           quot_options = list(),
           comp_fac = c(),
           dilution_options = list(),
           log_base = 2,
           do_impute = T,
           knn_options = list()) {
    stopifnot("SummarizedExperiment" %in% class(D))
    
    # function to merge default and user input
    map_lists <- function(def, user) {
      # check if there are any illegal entries
      bad_args <- setdiff(names(user), names(def))
      if (length(bad_args) > 0) {
        stop(sprintf(
          "Invalid argument(s): %s",
          paste0(bad_args, collapse = ", ")
        ))
      }
      # now fill up entries in def
      def[names(user)] <- user
      # return
      def
    }
    
    boxplot_options_def = list(
      plottile = "Sample boxplot",
      legend = T,
      ylabel = "Metabolite concentrations",
      logged = F,
      ggadd = NULL
    )
    
    dilution_options_def = list(boxpl = T, ggadd = NULL)
    knn_options_def = list(method = "knn.obs.euc.sel",
                           K = 10,
                           verbose = F)
    #Batch function must have paratemeter "D" and "batches"
    # add batch_ref_samples
    if (!missing(batch_col)) {
      #add another if statement
      if(!is.null(batch_ref_samples)){
        D %<>%
          batch_fun(batches = batch_col, ref_samples = dplyr::enquo(batch_ref_samples))
      } else{
        D %<>%
          batch_fun(batches = batch_col)
      }
    }

    boxplot_options <- map_lists(boxplot_options_def, boxplot_options)
    boxplot_options$D <- D
    D <- do.call("mt_plots_sampleboxplot", boxplot_options)
    
    quot_options_def = list(
      vars = 1:ncol(D),
      na_err = F,
      ref_samples = NULL,
      met_max = 1
    )
    quot_options <- map_lists(quot_options_def, quot_options)
    quot_options$D <- D
    if (!is.null(quot_options$ref_samples)) {
      D <- do.call("mt_pre_norm_quot", quot_options)
    } else{ # if NULL, remove the ref_sample
      quot_options$ref_samples <- NULL
      D <- do.call("mt_pre_norm_quot", quot_options)
    }

    if (!missing(comp_fac)) {
      dilution_options <-
        map_lists(dilution_options_def, dilution_options)
      
      for (s in comp_fac) {
        dilution_options$D <- D
        dilution_options$comp <- s
        D <- do.call("mt_plots_qc_dilutionplot", dilution_options)
      }
    }

    # check if there is any correlation between normalization factors and outcomes (bad sign if so)
    boxplot_options$D <- D
    D <- do.call("mt_plots_sampleboxplot", boxplot_options)
    # logging
    # base = 2
    D %<>%  mt_pre_trans_log(base = log_base)
    
    # KNN imputation
    if (do_impute == T) {
      knn_options <- map_lists(knn_options_def, knn_options)
      knn_options$D <- D
      D <- do.call("mt_pre_impute_knn", knn_options)
      boxplot_options$D <- D
      D <- do.call("mt_plots_sampleboxplot", boxplot_options)
    }
    D
  }

