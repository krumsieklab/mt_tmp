# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# KNN imputation
# K. Do based on Simone Wahls script from Missing values paper

# to use winner of MV paper do: imputeKNN(dat, methods="knn.obs.euc.sel", K=10)

# 18.12.17    - KD,JK
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



#' KNN imputation
#'
#' Code from Kiki Do's and Simone Wahl's paper
#' https://www.ncbi.nlm.nih.gov/pubmed/30830398
#'
#' To reproduce exact setting of winner method from paper:
#' mti_imputeKNN(dat, methods="knn.obs.euc.sel", K=10)
#'
#' Note: Most parts of this code are currently not used, since we only use the above-mentioned parameter setting in MT.
#'
#' @param dat input data matrix
#' @param methods which variation of KNN imputation to perform
#' @param cor.var.sel correlation threshold for neighbor selection
#' @param K number of neighbors to use
#' @param verbose output status messages?
#'
#' @returns imputed data matrix
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @noRd

mti_imputeKNN <- function(dat,
                   methods = c("knn.vars","knn.obs","knn.obs.euc","knn.obs.euc.sel"),
                   cor.var.sel = 0.2,
                   K=5,
                   verbose=T) {


  datimp <- dat
  ## knn.variable
  if("knn.vars" %in% methods) {
    incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
    if(verbose)message(paste0("Number of imcomplete variables: ", length(incom.vars)))

      # DMwR::knnImputation kann nicht damit umgehen, wenn weniger als k verf?gbare (obs) vars
      # also manuell:
      D2 <- as.matrix(stats::dist(t(scale(dat)),upper=T,diag=T)) # eucl
      D2[D2==0] <- NA
      if(verbose)message("Proceeding var... ")
      for (j in incom.vars) {
        if(verbose)cat(paste(j," "))
        comobs <- complete.cases(dat[,j])
        Mean <- mean(dat[,j],na.rm=T)
        SD <- sd(dat[,j],na.rm=T)
        if(any(!is.na(D2[,j]))) {
          KNNvars <- order(D2[,j],na.last=NA)
          KNNvars <- KNNvars[sapply(KNNvars, function(jj) any(!is.na(dat[!comobs,jj])))]
          } else KNNvars <- NULL

          dattmp <- dat
          KNNvars_sel <- KNNvars[1:min(K,length(KNNvars))]
          if(any(!is.na(D2[,j])) & length(KNNvars)>=1) dattmp[!comobs,j]  <-
            sapply(1:length(which(!comobs)),function(co) {
              dattmp_all <- scale(dattmp)[co,KNNvars]
              dattmp_sel <- scale(dattmp)[co,KNNvars_sel]
              if(any(!is.na(dattmp_sel))) ret <- sum((dattmp_sel*SD+Mean)*exp(-D2[KNNvars_sel,j]),na.rm=T)/sum(exp(-D2[KNNvars_sel,j])[!is.na(dattmp_sel)],na.rm=T) else ret <- na.omit(dattmp_all)[1]*SD+Mean
              ret})
          if(any(is.na(dattmp[,j]))) {
            still.incom <-  !complete.cases(dattmp[,j])
            dattmp[still.incom,j] <- mean(dattmp[!still.incom,j],na.rm=T) }
          datimp[,j] <- dattmp[,j]
      }

  } else if("knn.obs" %in% methods) { ## knn.sample: KNN per sample using mahalanobis distance (since dimensions=variables do not have same scale)

      # yaImpute::yai verwendet nur ganz complete obs! --> manuell!
      # dasselbe gilt f?r mahalanobis und mahalanobis.dist[StatMatch], wobei letzteres wenigstens ganze Matrix macht
      #--> manuell
      Sx <- cov(dat, use = "p")
      if(any(is.na(Sx))) {
        submatrix <- dat[-c(which(colSums(is.na(Sx))>0))]
        Sx <- cov(submatrix, use = "p")}
      incom.obs <- which(apply(dat,1,function(x) any(is.na(x))))
      if(verbose)message(paste0("Number of imcomplete observations: ", length(incom.obs)))
      if(verbose)message("Proceeding obs... ")
      for (i in incom.obs){
        if(verbose)cat(paste(i," "))
        comvars <-  complete.cases(as.numeric(dat[i,]))
        datpart <- dat[,comvars]
        obs_i <- as.numeric(datpart[i,])
        Sx_i <- Sx[colnames(datpart),colnames(datpart)]
        D2 <- sapply(1:nrow(dat),function(j) tryCatch(stats::mahalanobis(obs_i[!is.na(datpart[j,])],na.omit(as.numeric(datpart[j,])),Sx_i[!is.na(datpart[j,]),!is.na(datpart[j,])]),error=function(e) NA ))
        D2[D2==0] <- NA
        if(any(!is.na(D2))) {
          KNNids <- order(D2,na.last=NA)
          KNNids <- KNNids[sapply(KNNids,function(j) any(!is.na(dat[j,!comvars])))]
        }  else KNNids <- NULL


          dattmp <-  dat
          if(!is.null(KNNids)) KNNids_sel <- KNNids[1:min(K,length(KNNids))]
          if(any(!is.na(D2)) & length(KNNids)>=1) dattmp[i,!comvars] <-
            sapply(1:length(which(!comvars)), function(co) {
              dattmp_all <- dattmp[KNNids,co]
              dattmp_sel <- dattmp[KNNids_sel,co]
              if(any(!is.na(dattmp_sel))) ret <- sum(dattmp_sel*exp(-D2[KNNids_sel]),na.rm=T)/sum(exp(-D2[KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) else ret <- na.omit(dattmp_all)[1]
              ret})

          if(any(is.na(dattmp[i,]))) {        # weil z.B. D2 NAs enth?lt, ODER weil keine gemeinsam obs vars unter den KNN
            still.incom <-  !complete.cases(as.numeric(dattmp[i,]))
            dattmp[i,still.incom] <- apply(dattmp[,still.incom,drop=F],2,mean,na.rm=T)}   # wenn das nicht geht, mean ?ber alle vars

          datimp[i,] <- dattmp[i,]
          }

    } else if("knn.obs.euc" %in% methods) { ## knn.sample.euc: euclidean distance

      incom.obs <- which(apply(dat,1,function(x) any(is.na(x))))
      if(verbose)message(paste0("Number of imcomplete observations: ", length(incom.obs)))
      if(verbose)message("Proceeding obs... ")
      D2 <- as.matrix(stats::dist(scale(dat),upper=T,diag=T)) # eucl
      D2[D2==0] <- NA
      for (i in incom.obs){
        if(verbose)cat(paste(i," "))
        comvars <-  complete.cases(as.numeric(dat[i,]))
        if(any(!is.na(D2[i,]))) {
          KNNids <- order(D2[i,],na.last=NA)
          KNNids <- KNNids[sapply(KNNids,function(j) any(!is.na(dat[j,!comvars])))]
        }  else KNNids <- KNNids_naomit <- NULL


          dattmp <-  dat
          if(!is.null(KNNids)) KNNids_sel <- KNNids[1:min(K,length(KNNids))]
          if(any(!is.na(D2[i,])) & length(KNNids)>=1) dattmp[i,!comvars] <-
            sapply(1:length(which(!comvars)), function(co) {
              dattmp_all <- dattmp[KNNids,co]
              dattmp_sel <- dattmp[KNNids_sel,co]
              if(any(!is.na(dattmp_sel))) ret <- sum(dattmp_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) else ret <- na.omit(dattmp_all)[1]
              ret})

          if(any(is.na(dattmp[i,]))) {        # weil z.B. D2 NAs enth?lt, ODER weil keine gemeinsam obs vars unter den KNN
            still.incom <-  !complete.cases(as.numeric(dattmp[i,]))
            dattmp[i,still.incom] <- apply(dattmp[,still.incom,drop=F],2,mean,na.rm=T)}   # wenn das nicht geht, mean ?ber alle vars

          datimp[i,] <- dattmp[i,]
          }


    } else if("knn.obs.euc.sel" %in% methods) { ## knn.sample.euc.sel KNN per sample extended to cor cutoff (Tutz)  - use only variables with cor>0.2 for distance computation

      incom.obs <- which(apply(dat,1,function(x) any(is.na(x))))
      incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
      if(verbose)message(paste0("Number of imcomplete observations: ", length(incom.obs)))
      if(verbose)message("Proceeding obs... ")
      Cor <- cor(dat,use="p")
      D2list <- lapply(incom.vars, function(j) {
        varsel <- which(abs(Cor[j,])>cor.var.sel)    # ist j selbst dabei, ist aber ok
        if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
        if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
        D2 <- as.matrix(stats::dist(scale(dat[,varsel])),upper=T,diag=T)
        if(any(is.na(D2))) {
          D2a <- as.matrix(stats::dist(scale(dat)),upper=T,diag=T)*sqrt(length(varsel)/ncol(dat))
          D2[is.na(D2)] <- D2a[is.na(D2)] }
        diag(D2) <- NA
        D2})
      names(D2list) <- incom.vars
      for (i in incom.obs){
        if(verbose)cat(paste(i," "))
        comvars <-  complete.cases(as.numeric(dat[i,]))
        dattmp <-  dat
        for (j in which(!comvars)) {
          D2 <- D2list[[as.character(j)]]
          if(any(!is.na(D2[i,]))) {
            KNNids <- order(D2[i,],na.last=NA)
            KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(dat[ii,j])))]
          }  else KNNids  <- NULL

          # JK, 10/29/19
          # KNNids cannot actually be NULL, this happens if there's a row with all Infs or NA
          if (is.null(KNNids)) stop("Imputation method could not calculate correlations for some samples. Dataset probably contains samples that are all zero or NA and were then logged.")


            if(!is.null(KNNids)) KNNids_sel <- dplyr::intersect(KNNids[1:min(K,length(KNNids))],KNNids_naomit)
            if(length(KNNids_sel)<1) KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))] else
              if(length(which(sapply(KNNids_sel,function(ii) !is.na(dat[ii,j])))) < floor(K/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))]
            if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
              dattmp_sel <- dattmp[KNNids_sel,j]
              dattmp[i,j] <- sum(dattmp_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) }
        }
        datimp[i,] <- dattmp[i,]
        }


      datimp <- apply(datimp,2, function(x) {
          if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
          x})

    }


  datimp
  }
