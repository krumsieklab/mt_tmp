# QUESTIONS: DELETE ONCE ANSWERED
# 1. Is there a reason not to change sample_id to id_col?
# 2. Optional used a lot. When should optional be used?

# FIXES:
# 1. Add ability to accept SE object.
# 2. Change met_qc -> met_qc_sheet.
# 3. Change sample_qc -> sample_qc_sheet.
# 4. Change excel_file -> file.

#' Load Nightingale-format data
#'
#' Loads data from a Nightingale format Excel file.
#'
#' @param excel_file input excel file
#' @param format_type optional, type of nightingale format, single_sheet, multiple_sheets
#' @param data_sheet optional if format_type is multiple sheets, name of sheet with data
#' @param met_sheet optional, name of sheet with biomarker information
#' @param met_qc optional, name of sheet with biomarker qc tags
#' @param sample_qc optional, name of sheet with sample qc tags
#' @param sample_id optional, name of the column with sample id
#'
#' @return Produces an initial SummarizedExperiment, with assay, colData, rowData, and metadata with first entry.
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_files_load_nightingale(excel_file = codes.makepath("Mt/sampledata_night.xlsx")) %>%
#'   ...}
#'
#' @author RB
#'
#' @export
mt_load_nightingale <- function(excel_file,
                                data_sheet,
                                sample_id,
                                format_type = 'multiple_sheets',
                                met_sheet = 'Biomarker annotations',
                                met_qc = 'Tags per biomarker',
                                sample_qc = 'Quality control tags and notes') {
  if(format_type=='multiple_sheets'){
    if(missing(data_sheet)){data_sheet <- 'Results'}
    if(missing(sample_id)){sample_id <- 'Sample id'}
    D <- load_multiple_sheet_format(excel_file=excel_file,
                                    data_sheet=data_sheet, met_sheet=met_sheet,
                                    met_qc=met_qc, sample_qc=sample_qc, sample_id=sample_id)
  } else if (format_type=='single_sheet') {
    if(missing(data_sheet)){stop('data_sheet should be provided for single_sheet format!')}
    if(missing(sample_id)){sample_id <- 'sampleid'}
    D <- load_single_sheet_format(excel_file=excel_file,
                                    data_sheet=data_sheet, sample_id=sample_id)

  } else {stop(sprintf('Unknown format type, %s!', format_type))}

  # add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Nightingale file: %s, sheet: %s", basename(excel_file), data_sheet)
    )

  # return
  D

}

get_data <- function(mat, met_info, sample_id, sample_qc=F){
  # find last metabolite column
  imetlast <- max(which(apply(is.na(mat),2,sum)<dim(mat)[1]))
  # find last sample row
  isamplast <- max(which(apply(is.na(mat),1,sum)<dim(mat)[2]))
  # find row with table header
  tab_header <- min(which(!is.na(mat[, 1]))) + 1
  # find the first row of the data
  data_start <- min(which(!is.na(mat[(tab_header+1):isamplast, 1]))) + (tab_header + 1)
  # subset to data rows and columns
  if(sample_qc){
    col_ids <- c(which(mat[tab_header, ] %in%sample_id): imetlast)
    col_names <- unlist(c(sample_id, mat[(tab_header-1), col_ids[-1]]))
    mat <- mat[data_start:isamplast, col_ids]
    names(mat) <- col_names
  } else{
    # find the names of the table
    col_names <- c(sample_id, met_info[['Excel_column_name']])
    mat <- mat[data_start:isamplast, which(mat[tab_header, ] %in%col_names)]
    names(mat) <- col_names
  }
  return(mat)
}

load_multiple_sheet_format <- function(excel_file,
                                       data_sheet, met_sheet, met_qc, sample_qc, sample_id){
  # using readxl package:
  raw <- readxl::read_excel(path=excel_file, sheet=data_sheet, col_names = F)
  met_info <- readxl::read_excel(path=excel_file, sheet=met_sheet, col_names = T)
  qc_met <- readxl::read_excel(path=excel_file, sheet=met_qc, col_names = F)
  qc_sample <- readxl::read_excel(path=excel_file, sheet=sample_qc, col_names = F)
  # convert any spaces in the colnames to underscores
  colnames(met_info) <- gsub(" ", "_", colnames(met_info))

  result <- list()

  raw_data <- get_data(mat=raw, met_info, sample_id)
  qc_met <- get_data(mat=qc_met, met_info, sample_id) %>% t() %>% data.frame()
  names(qc_met) <- unlist(qc_met[1, ])
  qc_met <- qc_met[-1, ] %>% mutate('Excel_column_name' = rownames(qc_met[-1, ]))
  qc_sample <- get_data(mat=qc_sample, met_info, sample_id, sample_qc = T)
  # add sample information
  result$sampleinfo <- data.frame(raw_data %>% select(sample_id), stringsAsFactors = F, check.names = F) %>%
    left_join(qc_sample, by=sample_id)
  # add metabolite information
  result$metinfo <- data.frame(met_info, check.names = F) %>%
    left_join(qc_met, by='Excel_column_name')
  # add data
  raw_data <- raw_data %>% select(-sample_id) %>%
    mutate_all(as.matrix) %>% mutate_all(as.numeric)
  result$data <- data.frame(raw_data)

  # set info flags
  result$info$file <- excel_file
  result$info$raw_sheet <- raw_sheet
  result$info$met_sheet <- met_sheet
  result$info$met_qc <- met_qc
  result$info$sample_qc <- sample_qc

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
  # return SummarizedExperiment
  return(D)
}
load_single_sheet_format <- function (excel_file=excel_file,
                         data_sheet=data_sheet, sample_id=sample_id){
  # using readxl package:
  raw <- readxl::read_excel(path=excel_file, sheet=data_sheet, col_names = F)
  # find last metabolite column
  imetlast <- max(which(apply(is.na(raw),2,sum)<dim(raw)[1]))
  # find last sample row
  isamplast <- max(which(apply(is.na(raw),1,sum)<dim(raw)[2]))
  # find row with table header
  tab_header <- min(which(!is.na(raw[, 1]))) +1
  # find the first row of the data
  data_start <- 1 + (min(which(!is.na(raw[(tab_header+1):isamplast, 1]))) + (tab_header + 1))
  # subset to data rows and columns
  col_ids <- c(which(raw[tab_header, ] %in%sample_id): imetlast)
  col_names <- unlist(c(sample_id, raw[tab_header, col_ids[-1]]))
  raw_data <- raw[data_start:isamplast, col_ids]
  names(raw_data) <- col_names
  result <- list()
  # add sample information
  result$sampleinfo <- raw[data_start:isamplast, 1:col_ids[1]]
  names(result$sampleinfo) <- c(unlist(raw[(tab_header-1), 1:(col_ids[1]-1)]), sample_id)
  result$sampleinfo %<>% select(sample_id, everything())
  # add metabolite information
  result$metinfo <- data.frame(name=unlist(raw[tab_header, col_ids[-1]]),
                               fullname=unlist(raw[(tab_header-1), col_ids[-1]]),
                               check.names = F)
  # add data
  raw_data <- raw_data %>% select(-sample_id) %>%
    mutate_all(as.matrix) %>% mutate_all(as.numeric)
  result$data <- data.frame(raw_data)

  # set info flags
  result$info$file <- excel_file
  # fix variable names
  colnames(result$data) <- result$metinfo$name
  # generate summarized experiment
  D <- SummarizedExperiment(assay    = t(result$data),
                            colData  = result$sampleinfo,
                            rowData  = result$metinfo,
                            metadata = list(sessionInfo=utils::sessionInfo(), parseInfo=result$info))

  # ensure colnames and rownames exist
  if (is.null(colnames(D))) colnames(D) <- 1:ncol(D)
  if (is.null(rownames(D))) rownames(D) <- result$metinfo$name
  # return SummarizedExperiment
  return(D)
}
