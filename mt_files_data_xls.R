# MetaboTools
#
# Read data matrix from Excel sheet.
# That is, this will read the assay(), but not rowData(), or colData()
#
# last update: 2018-11-04
# authors: JK

require(readxl)
require(glue)
require(magrittr)
require(SummarizedExperiment)

mt_files_data_xls <- function(
  file,             # Excel file
  sheet,            # sheet name or number
  samplesInRows=T,  # read samples as rows (T) or as columns (F)
  ID                # if samplesInRows==T -> name of sample ID column... must be exactly one
  # if samplesInRows==F -> name of metabolite name column... must be exactly one
) {
  
  # load excel sheet
  df <- as.data.frame(read_excel(path=file,sheet=sheet,col_names=T))
  
  # ensure that sample ID column exists
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
  
  # construct SummarizedExperiment
  cd <- data.frame(as.character(colnames(assay)))
  colnames(cd)[1] <- ID
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

