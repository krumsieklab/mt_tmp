#' Remove redundant pathway annotations.
#'
#' Remove identical pathways from an existing SummarizedExperiment
#' data structure that has a column of pathway annotations.
#'
#' @param D \code{SummarizedExperiment} input
#' @param met_ID Column containing metabolite IDs
#' @param pw_ID Column containing pathways IDs
#'
#' @return rowData: redundant pathway annotation from SE pw_ID column will be filtered.
#'
#' @examples
#' # first annotate metabolites using smp_db and then remove redundant pathways
#' \dontrun{... %>%
#'   mt_anno_pathways_HMDB(in_col = "HMDb_ID", out_col = "smp_db",
#'   pwdb_name = "SMP", db_dir = codes.makepath("snippets/packages/metabotools_external/hmdb")) %>%
#'   mt_anno_pathways_remove_redundant(met_ID = "HMDb_ID", pw_ID = "smp_db") %>%
#' ...}
#'
#' @author Parviz Gomari
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_anno_pathways_remove_redundant <- function(
  D,                  # SummarizedExperiment input
  met_ID,         # Column containing metabolite IDs
  pw_ID              # Column containing pathways IDs
) {

  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if(!met_ID %in% names(rowData(D)))
    stop(glue::glue("met_ID is invalid. please enter an existing column name."))

  if(!pw_ID %in% names(rowData(D)))
    stop(glue::glue("pw_ID is invalid. please enter an existing column name."))

  row_data <-
    D %>%
    rowData() %>%
    as.data.frame() %>%
    dplyr::select(met_ID = !!met_ID,
                  pw_ID = !!pw_ID)

  mti_logstatus(glue::glue("creating grouping indices for the pathway IDs in {pw_ID} "))
  row_data_indexed <-
    row_data %>%
    dplyr::filter(!is.na(met_ID),
           pw_ID != "NULL") %>%
    tidyr::unnest(pw_ID) %>%
    dplyr::group_by(pw_ID) %>%
    dplyr::arrange(met_ID) %>%
    dplyr::mutate(met_IDs = stringr::str_c(met_ID, collapse = "|")) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(met_ID,
              pw_ID,
              pw_idx = dplyr::group_indices(., met_IDs))

  pw_ID_replacement <- row_data_indexed %>%
    dplyr::group_by(pw_idx) %>%
    dplyr::arrange(pw_ID) %>%
    dplyr::mutate(tmp_ID = dplyr::first(pw_ID),
           all_IDs = stringr::str_c(pw_ID %>% unique(), collapse = "|")) %>%
    dplyr::ungroup() %>%
    dplyr::select(met_ID,
                  ID = tmp_ID,
                  all_IDs) %>%
    dplyr::distinct() %>%
    dplyr::group_by(met_ID) %>%
    dplyr::summarise(ID = stringr::str_c(ID, collapse = ", ")) %>%
    dplyr::mutate(ID = stringr::str_split(ID, ", ")) %>%
    dplyr::right_join(row_data, by = "met_ID") %>%
    .$ID

  pwdb_summary_replacement <-
    dplyr::inner_join(metadata(D)$pathways[[pw_ID]],
               dplyr::select(row_data_indexed, -met_ID) %>%
                 dplyr::distinct(),
               by = c("ID" = "pw_ID")) %>%
    dplyr::group_by(pw_idx) %>%
    dplyr::arrange(ID) %>%
    dplyr::mutate(tmp_ID = dplyr::first(ID),
           all_IDs = stringr::str_c(ID, collapse = "|"),
           tmp_name = dplyr::first(pathway_name),
           all_pathway_names = stringr::str_c(pathway_name, collapse = "|")) %>%
    dplyr::filter(ID == tmp_ID) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(tmp_ID, tmp_name, pw_idx))

  # replace rowData and metadata of D
  rowData(D)[[pw_ID]] <- pw_ID_replacement
  metadata(D)$pathways[[pw_ID]] <- pwdb_summary_replacement

  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('removed redundant pathway annotations using the %s column', pw_ID)
    )

  D
}
