#' Load UC Davis-format data.
#'
#' Loads data from a UC Davis-format Excel file.
#'
#' @param D \code{SummarizedExperiment} input (missing if first step in pipeline)
#' @param file input Excel file in UCD format
#' @param sheet name or number of sheet
#' @param zero_to_NA replace zeros by NAs? (default: T)
#' @param gen_valid_varnames enforce valid R variable names?
#'
#' @return Produces an initial SummarizedExperiment, with assay, colData, rowData, and metadata with first entry
#'
#' @examples
#' \dontrun{D <-
#'   # load data
#'   mt_load_files_UCD(codes.makepath("Mt/sampledata.xlsx"), "OrigScale") %>%
#'   ...}
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_load_files_UCD <- function(
  D,
  file,                 # input xls file
  sheet,                # sheet name or number to read
  zero_to_NA=T,         # set zeros to NAs?
  gen_valid_varnames=F  # enforce valid variable names?
) {

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

  # using readxl package:
  raw = readxl::read_excel(path=file, sheet=sheet, col_names = F)

  #### find dimensions of sheet
  # find metabolite header row and last metabolite row
  imetheader = min(which(!is.na(raw[,1])))
  imetlast = max(which(apply(is.na(raw),1,sum)<dim(raw)[2]))
  # find sample header column and last sample row
  isampheader = min(which(!is.na(raw[1,])))
  isamplast = max(which(apply(is.na(raw),2,sum)<dim(raw)[1]))



  #### extract metabolite information
  # extract raw
  metinfo <- raw[imetheader:imetlast,1:(isampheader)]
  colnames(metinfo) <- metinfo[1,]
  metinfo <- metinfo[-1,]
  # find a metabolite name
  lookfor <- c("BinBase name","Metabolite name","Annotation")
  inds <- match(lookfor, colnames(metinfo))
  if (sum(!is.na(inds))==0) stop(sprintf("Could not find any metabolite name column. Was looking for %s", paste0(lookfor, collapse = ", ")))
  # pick first annotation
  ind <- inds[!is.na(inds)][1]
  metinfo$name <- metinfo[[ind]]
  # fix names
  if (gen_valid_varnames) colnames(metinfo) <- make.names(colnames(metinfo))
  # store
  result$metinfo <- metinfo


  #### extract sample information
  # extract raw
  sampleinfo <- data.frame(t(raw[1:(imetheader-1),(isampheader):isamplast]), stringsAsFactors = F)
  colnames(sampleinfo) <- sampleinfo[1,]
  if (gen_valid_varnames) colnames(sampleinfo) <- make.names(colnames(sampleinfo))
  sampleinfo <- sampleinfo[-1,]
  # store
  result$sampleinfo <- sampleinfo

  #### extract data
  result$data <- t(raw[(imetheader+1):imetlast, (isampheader+1):isamplast])
  result$data <- as.data.frame(apply(result$data,2, as.numeric))
  # set zeros to NAs?
  if (zero_to_NA) {
    result$data[result$data==0] = NA
  }


  #### return SummarizedExperiment, make sure colData and rowData are data.frames
  D <- SummarizedExperiment(assay    = t(result$data),
                            colData  = result$sampleinfo %>% as.data.frame(),
                            rowData  = result$metinfo %>% as.data.frame(),
                            metadata = list(sessionInfo=utils::sessionInfo(), parseInfo=result$info))

  # add original metadata if exists
  if (!is.null(result$meta$results)) metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) metadata(D)$settings <- result$meta$settings

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded UCD file: '%s', sheet: %s", basename(file), sheet)
    )

  # return
  D

}


# test code from jan
if (F) { # never execute automatically
  file = '/Users/jak2043/work/codes/neo4j/data/metMetabolon/10090_mmu/NonTargeted_Test sample_results_HELM-13-14ML+ CDT_150414.xlsx'
  sheet = 'OrigScale'

  D = parseMetabolonFile(file,sheet)
}
