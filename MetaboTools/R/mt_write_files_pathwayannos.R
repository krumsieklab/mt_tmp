#' Write out all pathway annotations, possible redundant pathways, and metabolite annotations.
#'
#' Creates an Excel file that contains detailed metabolite-to-pathway mapping information. Writes out 4 different forms of the mapping information.
#' Also includes information about which pathways are redundant, i.e. contain the same set of metabolites.
#'
#' @param D \code{SummarizedExperiment} input
#' @param pwfield name of the pathway annotation field
#' @param file output filename to write to
#'
#' @return Does not change the SummarizedExperiment.
#'
#' @examples
#' \dontrun{... %>%
#' mt_add_pathways(in_col = "KEGG",
#'     out_col = "kegg_db",
#'     pw_species = "hsapiens",
#'     pw_name = "kegg") %>%
#' mt_anno_pathways_remove_redundant(met_ID_col = "KEGG", pw_col = "kegg_db") %>%
#' mt_write_files_pathwayannos(pwfield='kegg_db', file='pwannos.xlsx')}
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_write_files_pathwayannos <- function(D,
                                        pwfield,
                                        file) {

  wb <- openxlsx::createWorkbook()
  # helper function
  add_sheet <- function(df, sheetname) {
    openxlsx::addWorksheet(wb, sheetname)
    openxlsx::writeData(wb, sheetname, df, rowNames = F, colNames=T)
  }

  ### pathways with one annotation per row ("long")
  # just generate, don't export yet
  pws.long <- D %>%
    metadata() %>%
    .$pathways %>%
    .[[pwfield]] %>%
    tidyr::separate_rows(all_IDs, all_pathway_names, sep='\\|')

  pws.wide <- D %>%
    metadata() %>%
    .$pathways %>%
    .[[pwfield]]


  ### metabolite annotations
  mets.long <- D %>%
    # build list of pathway IDs seperated by |
    rowData() %>%
    tibble::as.tibble() %>%
    dplyr::select(dplyr::one_of(c("name",pwfield))) %>%
    dplyr::rename(pw=dplyr::one_of(pwfield)) %>%
    dplyr::mutate(pw_str= sapply(pw, paste, collapse='|')) %>%
    dplyr::select(name, pw_str) %>%
    # explode into one row per annotation
    tidyr::separate_rows(pw_str, sep='\\|') %>%
    # add original pathway names
    dplyr::left_join(pws.wide, by = c("pw_str" = "ID")) %>%
    dplyr::select(name, all_IDs, all_pathway_names) %>%
    # explode redundant ones
    tidyr::separate_rows(all_IDs, all_pathway_names, sep='\\|') %>%
    # clean up
    dplyr::rename_( .dots=stats::setNames(list('all_IDs'), pwfield)) %>%
    dplyr::rename(pathway_name = all_pathway_names)
  # write out sorted by metabolites and by pathway_name
  mets.long %>% dplyr::arrange(name) %>% add_sheet(sheetname = 'fullmap_by_met')
  mets.long %>% dplyr::arrange(pathway_name) %>% add_sheet(sheetname = 'fullmap_by_pw')

  ### generate wide version of metabolites
  mets.wide <- mets.long %>%
    dplyr::group_by_(pwfield, 'pathway_name') %>%
    dplyr::summarise(metabolites = stringr::str_c(name, collapse = "|")) %>%
    dplyr::select(dplyr::one_of(pwfield,'metabolites'))

  ### export long pathway list
  pws.long %>%
    dplyr::mutate(is.redundant=  ifelse(duplicated(pws.long$pathway_name), 'yes', 'no')) %>%
    dplyr::left_join(mets.wide, by=c("ID"=pwfield)) %>%
    add_sheet(sheetname = 'pathways_long')

  ### export pathway annotations with redundant ones in list form ("wide")
  D %>%
    metadata() %>%
    .$pathways %>%
    .[[pwfield]] %>%
    dplyr::left_join(mets.wide, by=c("ID"=pwfield)) %>%
    add_sheet(sheetname = 'pathways_compact')


  # write out Excel file
  openxlsx::saveWorkbook(wb, file=file, overwrite=TRUE)


  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Pathway annotations in '%s' exported to Excel file '%s'", pwfield, file)
    )

  # pass SummarizedExperiment back, so pipeline can keep running
  D

}
