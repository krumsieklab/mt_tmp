# MetaboTools
#
# covariate adjustment of metabolites
#
# last update: 2019-03-29
# authors: MB
#
#! returns corrected data scaled
#

mt_pre_confounding_correction = function(
  D,       # SummarizedExperiment input
  formulae  # sample annotation column that contains batch info
) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = assay(D)
  
  # covariates of samples 
  dfc = colData(D)
  frm = as.formula(paste0(c("y",as.character(formulae)), collapse = ""))
  
  # run correction
  lm.fits <- apply(X %>% t,2, function(y) lm(frm, data.frame(y=y, dfc))) #residuals(%>% scale %>% t
  
  # overall effect of covariates to each variable as lm pvalue 
  p.lms  <- sapply(lm.fits, function(fit){
    f <- summary(fit)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    log(unname(p))
  })
  
  Xh = do.call(rbind, lapply(lm.fits, function(x) as.vector(scale( residuals(x) ))))
  rm(lm.fits)
  
 if(!identical(dim(X), dim(Xh))) 
   stop("missing values in confounding variables not supported")
  
  # make sure colnames and rownames are kept 
  colnames(Xh) <- colnames(X)
  rownames(Xh) <- rownames(X)
    
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("covariates adjustment: '%s'", as.character(formulae)),
      output = list(
        pvals   = p.lms,
        formula = formulae
      )
    )
  
  # return
  assay(D) = Xh
  D
  
}
