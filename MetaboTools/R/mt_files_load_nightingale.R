#' Load Nightingale-format data.
#'
#' Loads data from a Nightingale format Excel file. 
#'
#' @param excel_file input Excel file
#' @param raw_sheet name of sheet with raw data
#' @param met_sheet name of sheet with biomarker information
#' @param sample_id optional, name of the column with sample id
#'
#' @return Produces an initial SummarizedExperiment, with assay, colData, rowData, and metadata with first entry
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_files_load_nightingale(codes.makepath("Mt/sampledata_night.xlsx"), "Results") %>%
#'   ...}
#'
#' @author RB
#'
#' @importFrom magrittr %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_files_load_nightingale <- function(
  excel_file,           # Nightingale xls file
  raw_sheet, # sheet name or number to read
  met_sheet, # sheet name of number to read
  sample_id ='Sample id'# optional name of the column with sample id, default is 'Sample id'
) {

  # using readxl package:
  raw = readxl::read_excel(path=excel_file, sheet=raw_sheet, col_names = F)
  met_info = readxl::read_excel(path=excel_file, sheet=met_sheet, col_names = T)
  # convert any spaces in the colnames to underscores
  colnames(met_info) <- gsub(" ", "_", colnames(met_info))
  result=list()

  # find row with table header
  tab_header = min(which(!is.na(raw[, 1]))) + 1 
  # find the names of the table
  col_names <- c(sample_id, met_info[['Excel_column_name']])
  # find last metabolite column
  imetlast = max(which(apply(is.na(raw),2,sum)<dim(raw)[1]))
  # find last sample row
  isamplast = max(which(apply(is.na(raw),1,sum)<dim(raw)[2]))
  # find the first row of the data
  data_start <- min(which(!is.na(raw[(tab_header+1):isamplast, 1]))) + (tab_header + 1)

  # subset to data rows and columns 
  raw_data <- raw[data_start:isamplast, which(raw[tab_header, ] %in%col_names)]
  names(raw_data) <- col_names
  # add sample information
  result$sampleinfo = data.frame(raw_data %>% select(sample_id), stringsAsFactors = F)
  # add metabolite information
  result$metinfo <- as.data.frame(met_info)
  # add data
  raw_data <- raw_data %>% select(-sample_id) %>% 
    mutate_all(as.matrix) %>% mutate_all(as.numeric)
  result$data <- data.frame(raw_data)

  # set info flags
  result$info$file = excel_file
  result$info$raw_sheet = raw_sheet
  result$info$met_sheet = met_sheet

  # return SummarizedExperiment
  # add display name
  result$metinfo$name   <- result$metinfo$Excel_column_name
  # fix variable names
  colnames(result$data) <- result$metinfo$CSV_column_name
  # generate summarized experiment
  D <- SummarizedExperiment(assay    = t(result$data),
                       colData  = result$sampleinfo,
                       rowData  = result$metinfo,
                       metadata = list(sessionInfo=utils::sessionInfo(), parseInfo=result$info))

  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- 1:ncol(D)
  if (is.null(rownames(D))) rownames(D) <- result$metinfo$Excel_column_name

  # add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Nightingale file: %s, sheet: %s, and sheet: %s", basename(excel_file), raw_sheet, met_sheet)
    )

  # return
  D

}
