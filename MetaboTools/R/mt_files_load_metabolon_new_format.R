#' Load NEW Metabolon-format data
#'
#' Loads data from a NEW Metabolon-format Excel file. Needs to be in the original "Client Data Table" new format that they deliver.
#'
#' @param file input Excel file
#' @param raw_sheet name of sheet with raw data
#' @param met_sheet name of sheet with metabolite info
#' @param sample_sheet name of sheet with sample info
#'
#' @return Produces an initial SummarizedExperiment, with assay, colData, rowData, and metadata with first entry
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_files_load_metabolon_new_format(codes.makepath("Mt/sampledata.xlsx"), 
#'   raw_sheet="Raw Data", met_sheet="Chemical Annotation", sample_sheet="Sample Meta Data") %>%
#'   ...}
#'
#' @author RB
#'
#' @importFrom magrittr %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_files_load_metabolon_new_format <- function(
  file,           # Metabolon xls file
  raw_sheet,          # sheet name or number to read for raw data
  met_sheet,          # sheet name or number to read for metabolite info
  sample_sheet         # sheet name or number to read for sample info
) {

  cols2discard <- c('PARENT_SAMPLE_NAME', 'SAMPLE_TYPE', 'CLIENT_IDENTIFIER',
                    'MATRIX_DESIGNATION')
  # using readxl package:
  raw = readxl::read_excel(path=file, sheet=raw_sheet, col_names = T)
  met_info = readxl::read_excel(path=file, sheet=met_sheet, col_names = T)
  sample_info = readxl::read_excel(path=file, sheet=sample_sheet, col_names = T)

  result=list()
  result$data <- raw %>% select(-all_of(cols2discard))
  result$metinfo <- met_info[match(names(result$data), met_info$CHEM_ID), ]
  result$sampleinfo <- sample_info[match(raw$PARENT_SAMPLE_NAME, sample_info$PARENT_SAMPLE_NAME), ]
  
  # as column names, use "CHEMICAL_NAME", if available
  if ("CHEMICAL_NAME" %in% colnames(result$metinfo)) {
    colnames(result$data) = result$metinfo$CHEMICAL_NAME
  } else {
    colnames(result$data) = c()
  }

  # as row names, use "PARENT_SAMPLE_NAME", if available
  if ("PARENT_SAMPLE_NAME" %in% colnames(result$sampleinfo)) {
    rownames(result$data) = result$sampleinfo$PARENT_SAMPLE_NAME
  } else {
    rownames(result$data) = c()
  }

  # set info flags
  result$info$file = file
  result$info$raw_sheet = raw_sheet
  result$info$met_sheet = met_sheet
  result$info$sample_sheet = sample_sheet

  # return SummarizedExperiment

  # add display name
  result$metinfo$name   <- result$metinfo$CHEMICAL_NAME
  # fix variable names
  colnames(result$data) <- result$metinfo$CHEMICAL_NAME %>% make.names()
  # generate summarized experiment
  D <- SummarizedExperiment(assay    = t(result$data),
                       colData  = result$sampleinfo,
                       rowData  = result$metinfo,
                       metadata = list(sessionInfo=utils::sessionInfo(), parseInfo=result$info))

  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- 1:ncol(D)
  if (is.null(rownames(D))) rownames(D) <- result$metinfo$CHEMICAL_NAME

  # add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Metabolon file: %s, sheets: %s, %s, %s", 
                       basename(file), raw_sheet, met_sheet, sample_sheet)
    )

  # return
  D

}
