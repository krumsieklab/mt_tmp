#' Confounding correction using stepwise AIC
#'
#' \code{mt_pre_confounding_correction_stepwise_aic} returns the corrected data
#'
#' each metabolite is corrected by lm with the confounding variables determined by stepwise aic
#' and residuals are returned in lieu of uncorrected metabolites expressions.
#'
#' @param D \code{SummarizedExperiment} input
#' @param cols_to_cor vector of column numbers from colData to correct for
#' @param n_cores number of cores to use in parallelization
#'
#' @return SE with corrected data
#' @return $output: returns data.frame with the metabolite, its covars and the fit pval and rsq values
#'
#' @examples
#'  \dontrun{#... %>% mt_pre_confounding_correction_stepwise_aic(cols_to_cor = c(1, 4, 5), n_cores = 10)
#'  }
#'
#'
#' @author Annalise Schweickart, RB (modified on 2020-08-22)
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_pre_confounding_correction_stepwise_aic <- function(D, cols_to_cor, n_cores = 1) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.numeric(cols_to_cor))
  
  # there should be atleast one covariate information per sample
  Y <- D %>% colData() %>% data.frame()
  # names of covariates
  col_names <- names(Y)[cols_to_cor]
  # which samples have atleast one non NA covariate information ?
  non_na <- apply(Y [, cols_to_cor], 1, FUN=function(x) 
    length(which(is.na(x)))!=length(x)) %>% which()
  # how many have no covariate information at all?
  rem <- nrow(Y) - length(non_na)
  if(rem > 0) warning(sprintf("%d samples with no covariate info were removed!", rem))
  
  # subset the input based on missing covariate info
  if(length(non_na) > 0) {
    D <- D [, non_na]
    X <- t(assay(D))
    # bind covariate columns with metabolites for modelling
    model_data <- cbind.data.frame(X, colData(D)[, cols_to_cor])
  } else stop("No samples with any covariate info!")
  
  # per metabolite loop 
  outlist <- parallel::mclapply(1:ncol(X), FUN=function(i) {
    # all covariates and one metabolite
    form <- paste(colnames(X)[i], "~", paste(col_names, collapse=" + "))
    # compute coefficients
    mod <- stats::lm(stats::as.formula(form), data = model_data)
    # stepwise selection of covariates that have an effect on metabolite
    stepmod <- MASS::stepAIC(mod, direction = "backward", trace=0, k = log(nrow(model_data)))
    # selected covariates
    selected_covars <- as.character(stats::formula(stepmod)[3])
    # if no covariates are selected... 
    if(selected_covars==1){
      met_vec <- X[, i] # ...return the original metabolite value
      selected_covars <- fit_pval <- fit_rsq <- NA
    } else { # otherwise
      # fetch the names of the covariates
      selected_covars <- unlist(strsplit(selected_covars, split="+", fixed=TRUE))
      # collapse separated by semicolon
      selected_covars <- paste(gsub(" ", "", selected_covars), collapse=";")
      # get the significance of the model
      fit_pval <- signif(stats::pf(summary(stepmod)$fstatistic[1], summary(stepmod)$fstatistic[2], summary(stepmod)$fstatistic[3], lower.tail = FALSE), 3)
      # rsq of the model
      fit_rsq <- signif(summary(stepmod)$r.squared, 3)
      # return the residual --> corrected metabolite values
      met_vec <- stepmod$residuals + stepmod$coefficients["(Intercept)"] 
    }
    # result vector
    res <- list(met_vec, list(colnames(X)[i], selected_covars, fit_pval, fit_rsq))
    return(res)
  }, mc.cores= n_cores)
  
  #Bind outputs into assay for summarized experiment
  covar_adjusted <- do.call(rbind, lapply(outlist, function(x) unlist(x[[1]])))
  covars_log <- do.call(rbind, lapply(outlist, function(x) unlist(x[2])))
  colnames(covars_log) <-  c("metabolite", "covariates", "model.rsq", "model.pvalue")
  colnames(covar_adjusted) <- colnames(assay(D))
  rownames(covar_adjusted)<- rownames(assay(D))
  assay(D) <- covar_adjusted
  
  # add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('Adjusted for covariate effects with stepwiseAIC'),
      output = covars_log
    )
  D
}
