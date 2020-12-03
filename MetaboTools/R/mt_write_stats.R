# QUESTIONS: DELETE ONCE ANSWERED
# 1. Why is =NULL used here instead of just using missing()?

#' Export statistical results from pipeline into an Excel file
#'
#' Writes out the statistics data frame with p-values, adjusted p-values, fold changes, test statistics etc. into an Excel sheet.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output filename to write to.
#' @param stat_names Names of one or more statistical comparison to be written out. If NULL, will export all.
#' @param sort_by_p Automatically sort by p-values? Default: False.
#' @param met_name OPTIONAL. Name of rowdata column to include in the output.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{# Write out all results
#' ... %>% mt_write_stats(file="results.xlsx") %>%
#' # Write out specific result]
#' ... %>% mt_write_stats(file="results.xlsx", stat_names="comp1") %>%}
#'
#' @author JK, RB
#'
#' @export
mt_write_stats <- function(D, file, stat_names=NULL, sort_by_p=F, met_name=NULL) {

  # verify input arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.character(file))

  # get all stats entries
  S <- D %>% MetaboTools:::mti_res_get_stats_entries()
  allcomps <- S %>% purrr::map("output") %>% purrr::map("name") %>% unlist()

  # restrict to one or output all?
  if (!is.null(stat_names)) {
    # find the ones to filter to
    m <- match(stat_names, allcomps)
    # throw error if any are not found
    if (any(is.na(m))) {
      stop(sprintf("Cannot find the following stats entries to export: %s", paste0(stat_names[is.na(m)], collapse = ", ")))
    }
    # subset
    S <- S[m]
  }

  # export all
  wb <- openxlsx::createWorkbook()
  for (i in 1:length(S)) {
    # add to workbook
    df <- S[[i]]$output$table
    # if met_name is provided and it matches a column in row data add it to the output
    # else just silently skip this step
    if(is.null(met_name)==F && is.na(match(met_name, names(rowData(D))))==F){
        df <- cbind.data.frame(df, met_name=unlist(data.frame(rowData(D))[met_name]))
    }
    # sort?
    if (sort_by_p) {
      df %<>% arrange(p.value)
    }
    # add direction?
    if ("groups" %in% names(S[[i]]$output) && length(S[[i]]$output$groups)==2) {
      # generate vector of directions for indexing
      inds <- as.numeric(df$estimate>0)+1
      # translate into names
      df$dir_highin <- S[[i]]$output$groups[inds]
    }
    # output
    name <- S[[i]]$output$name
    ws=openxlsx::addWorksheet(wb,sheetName=name)
    openxlsx::writeData(wb=wb, sheet=name, x=df)
  }
  openxlsx::saveWorkbook(wb, file, overwrite=T)

  # add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Exported sheets '%s' to Excel file '%s'",
                       S %>% purrr::map("output") %>% purrr::map("name") %>% unlist() %>% paste0(collapse = ', '), file)
    )

  # pass SummarizedExperiment back, so pipeline can keep running
  D
}
