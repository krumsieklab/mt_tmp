#' Confounding correction 
#' 
#' \code{mt_pre_confounding_correction} returns the corrected data 
#' 
#' each metabolits is corrected by lm with the confounding variables and residuals are returned 
#' in lieu of uncorrected metabolites expressions.  
#'  
#' @param D SummarizedExperiment object 
#' @param formulae formula uncluding column names from sample annotation, i.e, ~batch+age 
#' @param strata strata classes for stratified confounding correction, should be a column from colData(D)
#' @param scalea should data be scaled after correaction
#'
#' @return corrected data together bit lm.fit pvalue per metabolites to show confounding affect on them
#' @export
#'
#' @examples
#' 
#'  # not run
#'  #... %>% mt_pre_confounding_correction( formulae = ~batch + age, strata = "RUN_DAY" )
#' 
#' @keywords confounding, lm
#' 
#' @author mubu
#' 
mt_pre_confounding_correction <- function(
  D,         # SummarizedExperiment input
  formulae,  # sample annotation column that contains batch info
  strata = NULL, # strata for stratified correction 
  scalea = F
) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = assay(D)
  
  # covariates of samples 
  dfc = colData(D)
  frm = as.formula(paste0(c("y",as.character(formulae)), collapse = ""))
  if(!is.null(strata)) strata = dfc[,strata]
  else strata = rep(1, nrow(dfc))
  
  # residual function which supports NA
  fresid <- function(fit){
    resids = fit$residuals 
    naa = fit$na.action
    if(is.null(naa)) return(resids)
    
    resh <- rep(NA, length(resids)+length(naa))
    names(resh) = seq(length(resh))
    names(resh)[naa] = names(naa)
    names(resh)[-naa] = names(resids)
    resh[-naa] = resids
    resh
  }
  
  # function that returns fit p value
  fpv<- function(fit){
    f <- summary(fit)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    log(unname(p))
  }
  
  # run correction
  lm.fits <- apply(X %>% t,2, function(y){ 
    pocs = lapply(unique(sort(strata)) %>% {names(.)=. ; .}, 
                  function(i) lm(frm, data.frame(y=y, dfc)[strata==i,]))
    
    # indices to go back original indices
    rh = unlist(lapply(pocs, fresid))[order(order(strata))]
    # model pvalues
    p.lms = sapply(pocs, fpv)
    list(rh=rh, p.lms=p.lms)
  }) #residuals(%>% scale %>% t
  
  # # overall effect of covariates to each variable as lm pvalue 
  # p.lms  <- sapply(lm.fits, function(fit){
  #   f <- summary(fit)$fstatistic
  #   p <- pf(f[1],f[2],f[3],lower.tail=F)
  #   log(unname(p))
  # })
  
 Xh = do.call(rbind, lapply(lm.fits, `[[`,"rh"))
 p.lms = sapply(lm.fits, `[[`,"p.lms")
  
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

# MetaboTools
#
# covariate adjustment of metabolites
#
# last update: 2019-03-29
# authors: MB
#
#! returns corrected data scaled
#
# consider missing values 
# has to be colnames

