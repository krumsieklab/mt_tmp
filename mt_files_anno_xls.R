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
  IDanno,    # column that contains ID information for mapping
  IDdata=IDanno, # column in data to map by
  nomaperr=F # throw error (T) or warning (F) if something does not map
) {
  
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(annosfor %in% c("samples","metabolites"))) stop("annosfor must be either 'samples' or 'metabolites'")
  
  # load excel sheet
  df <- as.data.frame(read_excel(path=file,sheet=sheet,col_names=T))
  # ensure that sample ID column exists
  if (!(IDanno %in% colnames(df))) stop(glue("sample ID column '{IDanno}' does not exist in '{basename(file)}, sheet '{sheet}'"))
  if (any(is.na(df[[IDanno]]))) stop(glue("sample ID column '{IDanno}' contains empty cells, '{basename(file)}, sheet '{sheet}'"))
  rownames(df) <- df[[IDanno]]
  
  # slightly different behavior for samples or metabolites
  if (annosfor=="samples") {
    # ensure the IDdata column exists
    if (!(IDdata %in% colnames(colData(D)))) stop(glue("ID column '{IDdata}' does not exist in current sample annotations of SE"))
    # check that all samples are found in the colnames of the existing dataset
    m <- match(df[[IDanno]], colData(D)[[IDdata]])
    if (any(is.na(m))) { 
      msg <- sprintf("The following sample IDs could not be found in the existing data matrix: %s",paste0(df[[IDdata]][is.na(m)],collapse=","))
      if (nomaperr) stop(msg)
      # else warning(msg)
    }
    # make everything a string
    df[[IDanno]] %<>% as.character()
    colData(D)[[IDdata]] %<>% as.character()
    # merge data frames
    newdf <- dplyr::left_join(data.frame(colData(D)), df, by = setNames(IDanno,IDdata))
    newdf[[IDanno]] <- newdf[[IDdata]] # make sure anno column name also exists (if different from data column name)
    stopifnot(all.equal(newdf[[IDdata]],colData(D)[[IDdata]])) # to make sure nothing was mixed up
    rownames(newdf) <- colnames(D)
    colData(D) <- DataFrame(newdf)
    
  } else if (annosfor=="metabolites") {
    # ensure the IDdata column exists
    if (!(IDdata %in% colnames(rowData(D)))) stop(glue("ID column '{IDdata}' does not exist in current metabolite annotations of SE"))
    # check that all metabolites are found in the $name column of the existing dataset
    m <- match(df[[IDanno]], rowData(D)[[IDdata]])
    if (any(is.na(m))) { 
      msg <- sprintf("The following metabolite IDs could not be found in the existing data matrix: %s",paste0(df[[IDdata]][is.na(m)],collapse=","))
      if (nomaperr) stop(msg)
      # else warning(msg)
    } 
    # make everything a string
    df[[IDanno]] %<>% as.character()
    rowData(D)[[IDdata]] %<>% as.character()
    # merge data frames
    newdf <- dplyr::left_join(data.frame(rowData(D)), df, by = setNames(IDanno,IDdata))
    newdf[[IDanno]] <- newdf[[IDdata]] # make sure anno column name also exists (if different from data column name)
    stopifnot(all.equal(newdf[[IDdata]],rowData(D)[[IDdata]])) # to make sure nothing was mixed up
    rownames(newdf) <- rownames(D)
    rowData(D) <- DataFrame(newdf)
  
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


