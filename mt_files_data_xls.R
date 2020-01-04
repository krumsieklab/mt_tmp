require(readxl)
require(glue)
require(magrittr)
require(SummarizedExperiment)

#' Load data matrix from Excel file.
#' 
#' Loads numerical data matrix from Excel sheet.
#' Can handle files with samples in rows or columns.
#' 
#' Default name for entity in columns will be 'name".
#'
#' @param file input Excel file 
#' @param sheet name or number of sheet
#' @param samplesInRows  read samples as rows (T) or as columns (F). default: T (samples as rows)
#' @param ID if samplesInRows==T -> name of sample ID column... must be exactly one
#'           if samplesInRows==F -> name of metabolite name column... must be exactly one
#'
#' @return Produces an initial SummarizedExperiment, with assay, colData, rowData, and metadata with first entry
#'
#' @examples
#' # Load data, two sheets with sample annotations, and one sheet with metabolite annotations from the same file
#' D <- 
#'   # load raw data
#'   mt_files_data_xls(file=file, sheet="data", samplesInRows=T, ID="SAMPLE_NAME") %>% 
#'   # sample annotations from metabolomics run
#'   mt_files_anno_xls(file=file, sheet="sampleinfo", annosfor="samples", IDanno = "SAMPLE_NAME") %>% 
#'   # sample annotations from clinical table
#'   mt_files_anno_xls(file=file, sheet="clinicals", annosfor="samples", IDanno="SAMPLE_NAME") %>% 
#'   # metabolite annotations`
#'   mt_files_anno_xls(file=file, sheet="metinfo", annosfor="metabolites", IDanno="BIOCHEMICAL", IDdata = "name") %>% 
#'   ...
#' 
#' @author JK
#' 
mt_files_data_xls <- function(file,
                              sheet,
                              samplesInRows = T,
                              ID) {
  
  
  # load excel sheet
  df <- as.data.frame(read_excel(path=file,sheet=sheet,col_names=T))
  
  # ensure that sample ID column exists
  if (missing(ID)) stop("No ID column provided")
  if (!(ID %in% colnames(df))) stop(glue("sample ID column '{ID}' does not exist in '{basename(file)}, sheet '{sheet}'"))
  # now convert to rownames
  df %<>% tibble::column_to_rownames(ID)  
  
  # construct assay
  if (samplesInRows) {
    # need to transpose
    assay = t(df)
  } else {
    assay = as.matrix(df)
  }
  
  # save original names as metabolite annotation, and ensure valid R name rownames
  metinfo = data.frame(name=rownames(assay))
  rownames(assay) %<>% make.names()
  
  # ensure that all columns are numeric
  rn <- rownames(assay)
  assay <- apply(assay,2,as.numeric)
  rownames(assay) <- rn
  
  # construct SummarizedExperiment
  cd <- data.frame(as.character(colnames(assay)))
  colnames(cd)[1] <- ifelse(samplesInRows,ID,'sample')
  D <- SummarizedExperiment(
    assay=assay,
    rowData=metinfo,
    colData=cd
  )
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue("loaded assay from Excel file '{basename(file)}, sheet '{sheet}'")
    )
  
  # return
  D
  
}

