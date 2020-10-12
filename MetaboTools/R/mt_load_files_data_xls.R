#' Load data matrix from Excel file.
#'
#' Loads numerical data matrix from Excel sheet.
#' Can handle files with samples in rows or columns.
#'
#' Default name for entity in columns will be 'name".
#'
#' @param D \code{SummarizedExperiment} input (missing if first step in pipeline)
#' @param file input Excel file
#' @param sheet name or number of sheet
#' @param samples_in_rows  read samples as rows (T) or as columns (F). default: T (samples as rows)
#' @param ID_col if samples_in_rows==T -> name of sample ID column... must be exactly one
#'           if samples_in_rows==F -> name of metabolite name column... must be exactly one
#' @param zero_to_NA replace zeros by NAs? (default: F)
#'
#' @return Produces an initial SummarizedExperiment, with assay, colData, rowData, and metadata with first entry
#'
#' @examples
#' # Load data, two sheets with sample annotations, and one sheet with metabolite annotations from the same file
#' \dontrun{D <-
#'   # load raw data
#'   mt_load_files_data_xls(file=file, sheet="data", samples_in_rows=T, ID_col="SAMPLE_NAME") %>%
#'   # sample annotations from metabolomics run
#'   mt_files_anno_xls(file=file, sheet="sampleinfo", annosfor="samples", IDanno = "SAMPLE_NAME") %>%
#'   # sample annotations from clinical table
#'   mt_files_anno_xls(file=file, sheet="clinicals", annosfor="samples", IDanno="SAMPLE_NAME") %>%
#'   # metabolite annotations`
#'   mt_files_anno_xls(file=file, sheet="metinfo", annosfor="metabolites", IDanno="BIOCHEMICAL", IDdata = "name") %>%
#'   ...}
#'
#' @author JK
#'
#' @importFrom magrittr %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_load_files_data_xls <- function(D,
                                  file,
                                  sheet,
                                  samples_in_rows = T,
                                  ID_col,
                                  zero_to_NA=F) {

  # initialize result list
  result=list()

  # validate arguments
  if (missing(file)) stop("file must be provided")
  if (missing(sheet)) stop("sheet must be provided")

  # save input information
  result$info$file <- file
  result$info$sheet <- sheet

  # get metadata from D if present
  if(!missing(D)){
    # validate SE
    if ("SummarizedExperiment" %in% class(D)) stop("D is not of class SummarizedExperiment")
    if (sum(assay(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  # load excel sheet
  df <- as.data.frame(readxl::read_excel(path=file,sheet=sheet,col_names=T))

  # ensure that sample ID_col column exists
  if (missing(ID_col)) stop("No ID column provided")
  if (!(ID_col %in% colnames(df))) stop(glue::glue("sample ID column '{ID_col}' does not exist in '{basename(file)}, sheet '{sheet}'"))
  # now convert to rownames
  df %<>% tibble::column_to_rownames(ID_col)

  # construct assay
  if (samples_in_rows) {
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

  # zeros to NAs?
  if (zero_to_NA) assay[assay==0] <- NA

  cd <- data.frame(as.character(colnames(assay)))
  colnames(cd)[1] <- ifelse(samples_in_rows,ID_col,'sample')

  # construct SummarizedExperiment
  D <- SummarizedExperiment(assay = assay,
                            rowData = metinfo,
                            colData = cd,
                            metadata = list(sessionInfo=utils::sessionInfo(), parseInfo=result$info))

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings


  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("loaded assay from Excel file '{basename(file)}, sheet '{sheet}'")
    )

  # return
  D

}

