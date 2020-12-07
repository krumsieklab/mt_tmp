#' Load data matrix from Excel file
#'
#' Loads numerical data matrix from Excel sheet. Can handle files with samples in rows or columns. Default name for entity in
#' columns will be "name".
#'
#' @param D \code{SummarizedExperiment} input. Missing if first step in pipeline.
#' @param file Name of input Excel file.
#' @param sheet Name or number of sheet.
#' @param samples_in_rows  Read samples as rows (T) or as columns (F). Default: T (samples as rows).
#' @param id_col If samples_in_rows==T -> name of sample ID column. If samples_in_rows==F -> name of metabolite name
#'    column. Must be exactly one.
#' @param zero_to_na Replace zeros by NAs? Default: F.
#'
#' @return If first step in pipeline, creates \code{SummarizedExperiment} object. Populates empty assay, colData, and rowData data frames.
#'
#' @examples
#' # Load data, two sheets with sample annotations, and one sheet with metabolite annotations from the same file
#' \dontrun{D <-
#'   # load raw data
#'   mt_load_data_xls(file=file, sheet="data", samples_in_rows=T, id_col="SAMPLE_NAME") %>%
#'   # sample annotations from metabolomics run
#'   mt_anno_xls(file=file, sheet="sampleinfo", anno_type="samples", anno_id = "SAMPLE_NAME") %>%
#'   # sample annotations from clinical table
#'   mt_anno_xls(file=file, sheet="clinicals", anno_type="samples", anno_id="SAMPLE_NAME") %>%
#'   # metabolite annotations`
#'   mt_anno_xls(file=file, sheet="metinfo", anno_type="metabolites", anno_id="BIOCHEMICAL", data_id = "name") %>%
#'   ...}
#'
#' @author JK
#'
#' @export
mt_load_data_xls <- function(D,
                             file,
                             sheet,
                             samples_in_rows = T,
                             id_col,
                             zero_to_na=F) {

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
    if ("SummarizedExperiment" %in% class(D) == F) stop("D is not of class SummarizedExperiment")
    if (length(assays(D))!=0) stop("Passed SummarizedExperiment assay must be empty!")

    # get metadata
    result$meta <- metadata(D)
  }

  # load excel sheet
  df <- as.data.frame(readxl::read_excel(path=file,sheet=sheet,col_names=T))

  # ensure that sample id_col column exists
  if (missing(id_col)) stop("No ID column provided")
  if (!(id_col %in% colnames(df))) stop(glue::glue("sample ID column '{id_col}' does not exist in '{basename(file)}, sheet '{sheet}'"))
  # now convert to rownames
  df %<>% tibble::column_to_rownames(id_col)

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
  if (zero_to_na) assay[assay==0] <- NA

  # save column names as sample annotation
  cd <- data.frame(as.character(colnames(assay)))
  colnames(cd)[1] <- ifelse(samples_in_rows,id_col,'sample')

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

