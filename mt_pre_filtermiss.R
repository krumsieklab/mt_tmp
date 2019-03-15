# MetaboTools
#
# Filter by missingness.
#
# Filters either samples or metabolites. Won't do both in one call.
#
# last update: 2019-3-15
# authors: JK, MB
#
# group option for threshold is added 
#
# todo: document output

# dependencies
library(tidyverse)
library(magrittr)
source(codes.makepath("MT/mt_internal_helpers.R"))

# function definition
mt_pre_filtermiss <- function(
  D,             # SummarizedExperiment input
  metMax=NA,     # maximum fraction of missing metabolites
  sampleMax=NA,   # maximum fraction of missing samples
  met_group = NA  # for each group of samples metMax be applied
) {
  
  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(!is.na(metMax) || !is.na(sampleMax))
  stopifnot(!(is.na(metMax) && is.na(sampleMax)))
  stopifnot(!(!is.na(metMax) && (metMax<0 || metMax>1)))
  stopifnot(!(!is.na(sampleMax) && (sampleMax<0 || sampleMax>1)))
  
  # perform filtering
  if (!is.na(metMax)) {
    NA.mat = is.na(assay(D))
    # metabolite
    if(is.na(met_group)){ # if group metmax
      na.stat = rowSums(NA.mat) 
      metKeep = na.stat/ncol(D) <= metMax
      D=D[metKeep,]
    }else{
      stopifnot(met_group %in% (colData(D) %>% colnames))
      xmet_group = colData(D)[,met_group]
      na.stat = xmet_group %>% unique %>% {. = as.character(.); names(.) = .; .} %>% 
        sapply(function(x) rowSums(NA.mat[, xmet_group == x])/sum(xmet_group == x))
      metKeep = rowSums( na.stat <= metMax ) == ncol(na.stat)
      D=D[metKeep,]
    }
  
    # add status information
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
      mti_generate_result(
        funargs = funargs,
        logtxt = sprintf('metabolites filtered, %d%%, %d of %d removed', metMax*100,sum(!metKeep),length(metKeep)),
        logshort = sprintf("filter met %d%%", metMax*100),
        output = list(kept=as.vector(metKeep), na.stat = na.stat[metKeep, ], na.mat = NA.mat[metKeep,])
      )
    
  } else {
    
    # sample
    sampleKeep = apply(is.na(assay(D)),2,sum)/nrow(D) <= sampleMax
    D=D[,sampleKeep]
    
    # add status information
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
      mti_generate_result(
        funargs = funargs,
        logtxt = sprintf('samples filtered, %d%%, %d of %d removed', sampleMax*100,sum(!sampleKeep),length(sampleKeep)),
        logshort = sprintf("filter samples %d%%", sampleMax*100),
        output = list(kept=as.vector(sampleKeep))
      )
    
  }
  
  # throw error if filtering caused empty dataset
  if (ncol(D)==0 || nrow(D)==0) stop("Filtering resulted in empty dataset.")
  
  # return
  D
}


