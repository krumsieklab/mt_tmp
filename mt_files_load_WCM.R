# MetaboTools
#
# Read WCM-format files.
# -> Expects metabolites as rows and an HMDB column.
#
# last update: 2018-12-10
# authors: JK


library(readxl)
library(stringr)
library(SummarizedExperiment)
library(tidyverse)

mt_files_load_WCM <- function(
  file,          # Metabolon xls file
  sheet,         # sheet name or number to read
  zeroToNA=T     # replace zeros by NAs?
) {
  
  # load file in raw format
  raw = as.data.frame(read_excel(path=file, sheet=sheet)) %>% column_to_rownames("compound")
  names = rownames(raw) 
  rownames(raw) = make.names(rownames(raw))
  # split off HMDB column if it exists
  if ("HMDB" %in% colnames(raw)) {
    metinfo = data.frame(name=names,HMDB=raw$HMDB)
    raw = raw %>% dplyr::select(-HMDB)
  } else {
    metinfo = data.frame(name=rownames(raw))
  }
  
  # zeros to NAs?
  if (zeroToNA) raw[raw==0] <- NA
  
  # generate summarized experiment
  D <- SummarizedExperiment(assay    = as.matrix(raw),
                            rowData  = metinfo,
                            colData  = data.frame(ID=colnames(raw)),
                            metadata = list(sessionInfo=sessionInfo()))
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded WCM file: %s, sheet: %s", basename(file), sheet)
    )
  
  # return
  D
  
}



