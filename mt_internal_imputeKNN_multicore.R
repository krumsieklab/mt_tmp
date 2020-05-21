# MT compatible function --------------------------------------------------

#' kNN multi-core internal method.
#'
#' 
#' Default settings are from the winning approach of Kiki's paper paper.
#' https://link.springer.com/article/10.1007%2Fs11306-018-1420-2
#' Specifically "knn.obs.euc.sel"
#' 
#' This script creates a tmp folder where distance matrices and imputed
#' values are stored
#
#' @param dat A matrix containing metabolites as columns and samples as rows
#' @param cor.var.sel correlation cutoff for metabolites, default: 0.2
#' @param K Number of nearest neighbors to consider, default is 10 (recommendation: do not touch)
#' @param mc_cores Number of cores to use for imputation, dafualt: 5
#' @param verbose T/F, whether to output intermediate steps, default: F
#'
#' @return assay: imputed data
#'
#' @examples
#' # in the context of a SE pipeline
#' ... %>% mt_pre_impute_knn() %>% ...    # standard call
#' 
#' @author Parviz Gomari
#' 
#' 


library(tidyverse)
library(parallel)

mt_internal_imputeKNN_multicore <- function(dat,
                                            cor.var.sel = 0.2,
                                            K=5,
                                            mc_cores = 5,
                                            verbose = T) {
  
  
  datimp <- dat
  # select row_index of missing observation
  incom_obs <- which(apply(dat,1,function(x) any(is.na(x)))) 
  
  # select column_index of missing observation
  incom_vars <- which(apply(dat,2,function(x) any(is.na(x)))) %>% names()
  
  
  if(verbose)message(paste0("Number of imcomplete observations: ", length(incom_obs)))
  
  # calculate correlation matrix using pairwise complete observations pairwise.complete.obs
  Cor <- cor(dat,use="p")
  
  
  dist_path <- "tmp/distance_matrices/"
  dir.create(dist_path, recursive = TRUE)
  
  # Precalculate values that will be reused in the mclapply below. Increases speed...
  scale_dat <- scale(dat)
  dist_scale_dat <- dist(scale_dat)
  
  # Get a list of sample distance matrices based on each incom_vars and varsel.
  # Basically what this is doing is that, if an incom.var is missing, find
  # metabolites with highest correlation, and based on these metabolites
  # calculate which samples are close.
  if(verbose)message("creating distance matrices")
  junk_collector <- mclapply(incom_vars, function(j) {
    
    # select variables (metabolites) with > cor.var.sel
    varsel <- which(abs(Cor[j,])>cor.var.sel)    # ist j selbst dabei, ist aber ok
    
    # If length varsel > 10, take index of the 11 highest ranking variables 
    # (metabolites) correlating with j
    if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
    
    # If length varsel < 5, take index of the 6 highest ranking variables 
    # (metabolites) correlating with j
    if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
    
    # calculate distance matrix of each sample based on scaled varsel data (column-wise)
    # Convert to a full matrix
    D2 <- as.matrix(dist(scale_dat[,varsel]),upper=T,diag=T) 
    
    # if D2 has missing values, replace them by distance matrix of whole dataframe,
    # scaled by length of varsel and ncol(dat)
    if(any(is.na(D2))) {
      D2a <- as.matrix(dist_scale_dat,upper=T,diag=T)*sqrt(length(varsel)/ncol(dat))
      D2[is.na(D2)] <- D2a[is.na(D2)] }
    diag(D2) <- NA
    saveRDS(D2, paste0(dist_path, j, ".rds"), compress = FALSE)
    NULL
  }, 
  mc.cleanup = TRUE, 
  mc.preschedule = TRUE, 
  mc.cores = mc_cores)
  
  
  
  impute_path <- "tmp/imputed_samples/"
  dir.create(impute_path, recursive = TRUE)
  
  if(verbose)message("calculating imputation values")
  junk_collector <- 
    mclapply(incom_obs, 
             function(i) {
               
               # for the observation/sample that has a missing metabolite,
               # select all metabolites that have values
               incomvars <-  dat[i, is.na(dat[i,])] %>% names()
               
               dattmp <-  dat
               
               # iterate over te missing metabolites
               for (j in incomvars) {
                 
                 # select sample distance matrix calculated for missing metabolite
                 D2_path <- paste0(dist_path, j, ".rds")
                 D2 <- readRDS(D2_path)     
                 
                 if(any(!is.na(D2[i,]))) {
                   
                   # order samples by decreasing distance
                   KNNids <- order(D2[i,],na.last=NA)
                   
                   # omit neighbors that also have the same metabolite j missing
                   KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(dat[ii,j])))]
                   
                 }  else KNNids  <- NULL
                 
                 # Select at most 10 nearest neighbors which overlap with KNNids_naomit
                 if(!is.null(KNNids)) KNNids_sel <- intersect(KNNids[1:min(K,length(KNNids))],KNNids_naomit)
                 
                 # If there really are no KNNids_sel, get at most floor(K/2) metbaolites from KNNids_naomit
                 # OR if there < floor(K/2) KNNids_sel metabolites, get at most floor(K/2) metbaolites 
                 # from KNNids_naomit
                 if(length(KNNids_sel)<1) KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))] else 
                   if(length(which(sapply(KNNids_sel,function(ii) !is.na(dat[ii,j])))) < floor(K/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))]
                 
                 # calculate the imputation value
                 if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
                   dattmp_sel <- dattmp[KNNids_sel,j]
                   dattmp[i,j] <- sum(dattmp_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) }
               }
               
               saveRDS(dattmp[i,], paste0(impute_path, i, ".rds"), compress = FALSE)
               
               NULL
             }, 
             mc.cleanup = TRUE, 
             mc.preschedule = TRUE, 
             mc.cores = mc_cores)
  
  
  
  
  # Fill in imputed values
  
  for (i in incom_obs){
    
    datimp[i,] <- readRDS(paste0(impute_path, i, ".rds"))
    
  }
  
  
  # if there are still values missing, just assign mean value
  datimp <- apply(datimp,2, function(x) {
    if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
    x})
  
  
  datimp
}






