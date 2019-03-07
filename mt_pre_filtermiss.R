# MetaboTools
#
# Filter by missingness.
#
# Filters either samples or metabolites. Won't do both in one call.
#
# last update: 2018-10-11
# authors: JK
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
  sampleMax=NA   # maximum fraction of missing samples
) {
  
  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(!is.na(metMax) || !is.na(sampleMax))
  stopifnot(!(is.na(metMax) && is.na(sampleMax)))
  stopifnot(!(!is.na(metMax) && (metMax<0 || metMax>1)))
  stopifnot(!(!is.na(sampleMax) && (sampleMax<0 || sampleMax>1)))
  
  # perform filtering
  if (!is.na(metMax)) {
    
    # metabolite
    metKeep = apply(is.na(assay(D)),1,sum)/ncol(D) <= metMax
    D=D[metKeep,]
    # add status information
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
      mti_generate_result(
        funargs = funargs,
        logtxt = sprintf('metabolites filtered, %d%%, %d of %d removed', metMax*100,sum(!metKeep),length(metKeep)),
        logshort = sprintf("filter met %d%%", metMax*100),
        output = list(kept=as.vector(metKeep))
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


