#' Load Nightingale-format data.
#'
#' Loads data from a Nightingale format Excel file. 
#'
#' @param excel_file input Excel file
#' @param raw_sheet name of sheet with raw data
#' @param met_sheet name of sheet with biomarker information
#' @param met_qc name of sheet with biomarker qc tags
#' @param sample_qc name of sheet with sample qc tags
#' @param sample_id optional, name of the column with sample id
#'
#' @return Produces an initial SummarizedExperiment, with assay, colData, rowData, and metadata with first entry
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_files_load_nightingale(excel_file = codes.makepath("Mt/sampledata_night.xlsx"), 
#'   raw_sheet = "Results", 
#'   met_sheet = "Biomarker annotations",
#'   met_qc = "Tags per biomarker",
#'   sample_qc = "Quality control tags and notes",
#'   sample_id = 'Sample id'# optional name of the column with sample id, default is 'Sample id'
#'   ) %>%
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
  raw_sheet = 'Results', # sheet name or number to read
  met_sheet = 'Biomarker annotations', # sheet name of number to read
  met_qc = 'Tags per biomarker', 
  sample_qc = 'Quality control tags and notes',
  sample_id = 'Sample id'# optional name of the column with sample id, default is 'Sample id'
) {

  # using readxl package:
  raw <- readxl::read_excel(path=excel_file, sheet=raw_sheet, col_names = F)
  met_info <- readxl::read_excel(path=excel_file, sheet=met_sheet, col_names = T)
  qc_met <- readxl::read_excel(path=excel_file, sheet=met_qc, col_names = F)
  qc_sample <- readxl::read_excel(path=excel_file, sheet=sample_qc, col_names = F)
  # convert any spaces in the colnames to underscores
  colnames(met_info) <- gsub(" ", "_", colnames(met_info))
  
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
