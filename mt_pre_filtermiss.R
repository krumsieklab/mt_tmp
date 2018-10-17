# MetaboTools
#
# Filter by missingness.
#
# Filters either samples or metabolites. Won't do both in one call.
#
# last update: 2018-10-11
# authors: JK
#

# dependencies
library(tidyverse)
library(magrittr)
source(codes.makepath("packages/metabotools/mt_internal_helpers.R"))

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
  # perform filtering
  if (!is.na(metMax)) {
    
    # metabolite
    metKeep = apply(is.na(assay(D)),1,sum)/ncol(D) <= metMax
    D=D[metKeep,]
    # add status information
    metadata(D)$preprocess %<>% 
      add_to_list(list(txt=sprintf('metabolites filtered, %d%%, %d of %d removed', metMax*100,sum(!metKeep),length(metKeep)),call=match.call(), op="filter",dir="met",cutoff=metMax,kept=as.vector(metKeep)))
  } else {
    
    # sample
    sampleKeep = apply(is.na(assay(D)),2,sum)/nrow(D) <= sampleMax
    D=D[,sampleKeep]
    # add status information
    metadata(D)$preprocess %<>% 
      add_to_list(list(txt=sprintf('samples filtered, %d%%, %d of %d removed', sampleMax*100,sum(!sampleKeep),length(sampleKeep)), call=match.call(), op="filter",dir="samples",cutoff=sampleMax,kept=as.vector(sampleKeep)))
    
  }
  # return
  D
}


