library(MASS)
#' Confounding correction using stepwise AIC
#' 
#' \code{mt_pre_confounding_correction_stepwise_aic} returns the corrected data 
#' 
#' each metabolite is corrected by lm with the confounding variables determined by stepwise aic 
#' and residuals are returned in lieu of uncorrected metabolites expressions.  
#'  
#' @param D SummarizedExperiment object 
#' @param id_col column name of summarized experiment with sample id
#' @param med_file name of file containing medication data
#' @param to_remove vector of medication names that should not be corrected for
#' @param mc.cores number of cores to use in parallelization
#'
#' @return SE with corrected data 
#' @return $output: returns data.frame with the metabolite, its covars and the fit pval and rsq values
#' 
#' @examples
#' 
#'  #... %>% mt_pre_confounding_correction_stepwise_aic(id_col = "RID", med_file = meds, to_remove = c( "Med.Anticholinesterases", "Med.NMDAAntag"), mc.cores = 10)
#' 
#' 
#' @author Annalise Schweickart
#' 

mt_pre_confounding_correction_stepwise_aic <- function(
  D,         # SummarizedExperiment input
  id_col,    # Column name of ID for medication/sample matching
  med_file,  # File name containing medication data
  to_remove = c(), # any columns in med_file that should not be included in the correction
  mc.cores = 1  # number of cores to use in parallelization
) {
  ## Load medication data, remove Alzheimer's related meds
  meds <- read.table(med_file, header = T, sep = "\t")
  meds <- meds[,!colnames(meds)%in%to_remove]
  
  X <- data.frame(cbind(colData(D)[[id_col]], t(assay(D))))
  colnames(X)<- c(id_col, colnames(t(assay(D))))
  model.data <- merge(X, meds, by=id_col)
  meds.log <- NULL
  
  #parallelization step
  outlist <- mclapply(colnames(model.data[2:ncol(X)]), function(met) single_met_covars(
    met,
    id_col,
    meds,
    model.data,
    X
  ), mc.cores= mc.cores)
  
  #Bind outputs into assay for summarized experiment
  med.adjusted <- do.call(rbind, lapply(outlist, function(x) unlist(x[1])))
  meds.log <- do.call(rbind, lapply(outlist, function(x) unlist(x[2])))
  colnames(meds.log) <-  c("metabolite", "medications", "model.rsq", "model.pvalue")
  colnames(med.adjusted) <- colnames(assay(D))
  rownames(med.adjusted)<- rownames(assay(D))
  assay(D) <- med.adjusted
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('Adjusted for medication effects'),
      output = meds.log
    )
  
  D
}

single_met_covars <- function(
  metabolite,
  id_col,
  meds,
  model.data,
  X
){
  ## Select medications via backward-selection
  form <- paste(metabolite, "~", paste(colnames(meds)[2:ncol(meds)], collapse="+"))
  mod <- lm(as.formula(form), data = model.data)
  stepmod <- stepAIC(mod, direction = "backward", trace=0, k = log(nrow(model.data)))
  selected.covars <- as.character(formula(stepmod)[3])
  
  ## if no medications are selected
  if(selected.covars=="1"){
    selected.covars <- NA
    fit.pval <- NA
    fit.rsq <- NA
    met.vec <- X[, metabolite]
  }else{
    selected.covars <- unlist(strsplit(selected.covars, split="+", fixed=TRUE))
    selected.covars <- paste(gsub(" ", "", selected.covars), collapse=";")
    fit.pval <- signif(pf(summary(stepmod)$fstatistic[1], summary(stepmod)$fstatistic[2], summary(stepmod)$fstatistic[3], lower.tail = FALSE), 3)
    fit.rsq <- signif(summary(stepmod)$r.squared, 3)
    met.vec <- stepmod$residuals[match(X[[id_col]], model.data[[id_col]])] + stepmod$coefficients["(Intercept)"]
  }
  list(met.vec, list(metabolite, selected.covars, fit.pval, fit.rsq))
}
