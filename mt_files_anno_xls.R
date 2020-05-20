library(readxl)
library(glue)
library(magrittr)
library(SummarizedExperiment)

#' Load annotations from Excel file.
#' 
#' Loads annotations and merges them into current SummarizedExperiment. 
#' Performs "left-joins", i.e. leaves the original SE unchanged and just adds information where it can be mapped.
#' Can load annotations for both metabolites (rowData) and samples (colData)
#'
#' @param D \code{SummarizedExperiment} input
#' @param file input Excel file
#' @param sheet name or number of sheet
#' @param annosfor "samples" or "metabolites"
#' @param IDanno column in annotation file that contains ID information for mapping
#' @param IDdata column in existing data for mapping
#' @param nomaperr # throw error (T) or warning (F) if something does not map. default: F
#' @param colnames.from # take column from new annotation file and overwite colnames (i.e. sample names) of SE (default: none)
#'
#' @return rowData or colData: new annotations added
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
mt_files_anno_xls <-
  function(D,
           file,
           sheet,
           annosfor,
           IDanno,
           IDdata = IDanno,
           nomaperr = F,
           colnames.from = NA) {
    
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
    newdf <- dplyr::left_join(data.frame(colData(D), check.names=F), df, by = setNames(IDanno,IDdata))
    newdf[[IDanno]] <- newdf[[IDdata]] # make sure anno column name also exists (if different from data column name)
    stopifnot(all.equal(newdf[[IDdata]],colData(D)[[IDdata]])) # to make sure nothing was mixed up
    rownames(newdf) <- colnames(D)
    colData(D) <- DataFrame(newdf)
    
    # colnames.from?
    if (!is.na(colnames.from)) {
      # copy field, set colnames
      cn = newdf[[colnames.from]]
      colnames(D) = cn
    }
    
    
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
    
    # colnames.from?
    if (!is.na(colnames.from)) {
      stop("colnames.from cannot be used for metabolite annotations")
    }
  
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


