#' Export statistical results from pipeline into Excel file.
#' 
#' Writes out the statistics data frame with p-values, adjusted p-values, fold changes, test statistics etc. into an Excel sheet.
#'
#' @param D \code{SummarizedExperiment} input
#' @param file output filename to write to
#' @param compnames Names of one or more statistical comparison to be written out. If NULL, will export all.
#' @param sort.by.p Automatically sort by p-values? (default: False)
#'
#' @return Does not change the SummarizedExperiment.
#'
#' @examples
#' \dontrun{# Write out all results
#' ... %>% mt_files_write_stats(file="results.xlsx") %>%
#' # Write out specific result]
#' ... %>% mt_files_write_stats(file="results.xlsx", compnames="comp1") %>%}
#' 
#' @author JK
#' 
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#' 
#' @export
mt_files_write_stats <- function(D, file, compnames=NULL, sort.by.p=F) {
  
  # verify input arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.character(file))
  
  # get all stats entries
  S <- D %>% mti_res_get_stats_entries()
  allcomps <- S %>% purrr::map("output") %>% purrr::map("name") %>% unlist()
  
  # restrict to one or output all?
  if (!is.null(compnames)) {
    # find the ones to filter to
    m <- match(compnames, allcomps)
    # throw error if any are not found
    if (any(is.na(m))) { 
      stop(sprintf("Cannot find the following stats entries to export: %s", paste0(compnames[is.na(m)], collapse = ", ")))
    }
    # subset
    S <- S[m]
  }
  
  # export all
  wb <- openxlsx::createWorkbook()
  for (i in 1:length(S)) {
    # add to workbook
    df <- S[[i]]$output$table
    # sort?
    if (sort.by.p) {
      df %<>% arrange(p.value)
    }
    # output
    name <- S[[i]]$output$name
    ws=openxlsx::addWorksheet(wb,sheetName=name)
    openxlsx::writeData(wb=wb, sheet=name, x=df)
  }
  openxlsx::saveWorkbook(wb, file, overwrite=T)
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Exported sheets '%s' to Excel file '%s'", 
                       S %>% purrr::map("output") %>% purrr::map("name") %>% unlist() %>% paste0(collapse = ', '), file)
    )
  
  
  # pass SummarizedExperiment back, so pipeline can keep running
  D
}