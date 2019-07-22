library(openxlsx)

#' Outputs assay, colData and rowData into an Excel file.
#'
#' @param D \code{SummarizedExperiment} input
#' @param file output filename to write to
#'
#' @return Does not change the SummarizedExperiment.
#'
#' @examples
#' %>% mt_files_write_xls(file = "out.xlsx") %>%
#' @author JK, BGP
#' 
mt_files_write_xls <- function(D, file) {
  
  # verify that input is a SummarizedExperiment
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.character(file))
  
  # write out
  wb <- createWorkbook()
  addWorksheet(wb,"assay")
  writeData(wb, "assay", assay(D), rowNames = T, colNames=T)
  addWorksheet(wb,"rowData")
  writeData(wb, "rowData", rowData(D) %>% as.data.frame(), rowNames = T)
  addWorksheet(wb,"colData")
  writeData(wb, "colData", colData(D) %>% as.data.frame())
  saveWorkbook (wb, file=file, overwrite=TRUE) 
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Data exported to Excel file '%s'", file)
    )
  
  
  # pass SummarizedExperiment back, so pipeline can keep running
  D
  
}
