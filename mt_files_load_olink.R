require(readxl)
require(stringr)
require(SummarizedExperiment)
require(tidyverse)
require(OlinkAnalyze)

#' Load Olink-format data.
#' 
#' Loads data from an Olink-format Excel file.  
#' Uses Olink's R code from https://github.com/Olink-Proteomics/OlinkRPackage/tree/master/OlinkAnalyze
#' To install Olink code use:
#' install.packages("devtools")
#' devtools::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze')
#' 
#' In case the Olink file is not in XLSX format, but CSV or TSV text format:
#' you may need to remove all Ctrl data columns
#' save file as xlsx (make sure to rename SEPT9 back from Excel date)
#' don't keep the Excel file open while running R - this throws a file access denied error
#'
#' @param filename filename of NPX file exported from NPX manger.
#'
#' @return Produces an initial SummarizedExperiment, with assay (note: 2^NPX), colData, rowData, and metadata with first entry
#'
#' @examples
#' D <- 
#'   # load data
#'   mt_files_load_olink(codes.makepath("mt/olink_sampledata.xlsx")) %>%
#'   ...
#' 
#' @author KS
#' 
mt_files_load_olink <- function(
  file           # Olink Metabolon xlsx file
) {

  result=list()
  result$info = "nothing yet"
  
  # read the Olink file  
  odata = OlinkAnalyze::read_NPX(file)

  # NPX values are ~ PCA cycles, hence on some kind of a log-scale
  # to be in-line with meta-tools, un-logscale the data using 2 as a basis
  odata$NPX = 2^odata$NPX
  
  # convert to wide format
  wdata = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "NPX")
  xdata = simplify2array(wdata[,-1])
  rownames(xdata) = wdata$SampleID
  
  # convert to SummarizedExperiment
  D = SummarizedExperiment(assay = list(exprs = t(xdata)),
                              colData = DataFrame(sample_id = rownames(xdata)),
                              rowData = DataFrame(feature_id = colnames(xdata),
                                                  LODdata = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "LOD")[1,-1] %>% unlist() %>% as.numeric() %>% unname(),
                                                  MissingFreq = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "MissingFreq")[1,-1] %>% unlist() %>% as.numeric() %>% unname(),
                                                  OlinkID = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "OlinkID")[1,-1] %>% unlist() %>% unname(),
                                                  UniProt = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "UniProt")[1,-1] %>% unlist() %>% unname(),
                                                  Assay = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "Assay")[1,-1] %>% unlist() %>% unname(),
                                                  Panel = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "Panel")[1,-1] %>% unlist() %>% unname(),
                                                  Panel_Version = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "Panel_Version")[1,-1] %>% unlist() %>% unname(),
                                                  PlateID = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "PlateID")[1,-1] %>% unlist() %>% unname()
                              ),
                              metadata = list(sessionInfo=sessionInfo(), parseInfo=result$info)
  )
  
  # do some sanity checks
  stopifnot (length(which(colnames(xdata) != rowData(D)$OlinkID)) == 0)
  MissingFreq2 = reshape2::dcast(odata,formula = " SampleID ~ OlinkID", value.var = "MissingFreq")[2,-1] %>% unlist() %>% as.numeric() %>% unname()
  stopifnot (length(which(MissingFreq2 != rowData(D)$MissingFreq)) == 0)
  
  # as column names, use "BIOCHEMICAL", if available
  rowData(D)$BIOCHEMICAL = rowData(D)$Assay
  
  # as row names, use "SAMPLE_NAME", if available 
  D$SAMPLE_NAME = D$sample_id

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Olink file: %s", basename(file))
    )
  
  # return
  D

}
