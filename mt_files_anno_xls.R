# MetaboTools
#
# Add annotations (either sample or metabolite) from Excel file.
# Will add to colData(), or rowData()s
#
# last update: 2018-11-05
# authors: JK

require(readxl)
require(glue)
require(magrittr)
require(SummarizedExperiment)

mt_files_anno_xls <- function(
  D,         # SummarizedExperiment input
  file,      # Excel file
  sheet,     # sheet name or number
  annosfor,  # "samples" or "metabolites"
  ID         # column that contains ID information for mapping
) {
  
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(annosfor %in% c("samples","metabolites"))) stop("annosfor must be either 'samples' or 'metabolites'")
  
  # load excel sheet
  df <- as.data.frame(read_excel(path=file,sheet=sheet,col_names=T))
  # ensure that sample ID column exists
  if (!(ID %in% colnames(df))) stop(glue("sample ID column '{ID}' does not exist in '{basename(file)}, sheet '{sheet}'"))
  if (any(is.na(df[[ID]]))) stop(glue("sample ID column '{ID}' contains empty cells, '{basename(file)}, sheet '{sheet}'"))
  rownames(df) <- df[[ID]]
  
  # slightly different behavior for samples or metabolites
  if (annosfor=="samples") {
    # check that all samples are found in the colnames of the existing dataset
    m <- match(df[[ID]], colnames(D))
    if (any(is.na(m))) stop(sprintf("The following sample IDs could not be found in the existing data matrix: %s",paste0(df[[ID]][is.na(m)],collapse=",")))
    # add data
    newdf <- DataFrame(merge(data.frame(colData(D)), df, by='row.names', all.x=T) %>% select(-Row.names))
    rownames(newdf) <- colnames(D)
    colData(D) <- newdf
    
  } else if (annosfor=="metabolites") {
    # check that all metabolites are found in the $name column of the existing dataset
    m <- match(df[[ID]], rowData(D)$name)
    if (any(is.na(m))) stop(sprintf("The following metabolite IDs could not be found in the existing data matrix: %s",paste0(df[[ID]][is.na(m)],collapse=",")))
    # add data
    newdf <- merge(data.frame(rowData(D)), df, by.x='name', by.y='row.names', all.x=T)  
    rownames(newdf) <- rownames(D)
    rowData(D) <- newdf
    
  } else
    stop("bug")
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue("loaded {annosfor} annotations from Excel file '{basename(file)}, sheet '{sheet}'")
    )
  
  # return
  D
  
  
}


