#' Remove redundant pathway annotations
#'
#' Remove identical pathways from an existing SummarizedExperiment
#' data structure that has a column of pathway annotations.
#'
#' @param D \code{SummarizedExperiment} input
#' @param met_id Column containing metabolite IDs
#' @param pw_id Column containing pathways IDs
#'
#' @return rowData: Redundant pathway annotation from SE pw_id column filtered.
#'
#' @examples
#' # first annotate metabolites using smp_db and then remove redundant pathways
#' \dontrun{... %>%
#'   mt_anno_pathways_HMDB(in_col = "HMDb_ID", out_col = "smp_db",
#'   pwdb_name = "SMP", db_dir = codes.makepath("snippets/packages/metabotools_external/hmdb")) %>%
#'   mt_anno_pathways_remove_redundant(met_id = "HMDb_ID", pw_id = "smp_db") %>%
#' ...}
#'
#' @author PG
#'
#' @export
mt_anno_pathways_remove_redundant <- function(D,
                                              met_id,
                                              pw_id ) {

  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if(!met_id %in% names(rowData(D)))
    stop(glue::glue("met_id is invalid. please enter an existing column name."))

  if(!pw_id %in% names(rowData(D)))
    stop(glue::glue("pw_id is invalid. please enter an existing column name."))

  row_data <-
    D %>%
    rowData() %>%
    as.data.frame() %>%
    dplyr::select(met_id = !!met_id,
                  pw_id = !!pw_id)

  mti_logstatus(glue::glue("creating grouping indices for the pathway IDs in {pw_id} "))
  row_data_indexed <-
    row_data %>%
    dplyr::filter(!is.na(met_id),
           pw_id != "NULL") %>%
    tidyr::unnest(pw_id) %>%
    dplyr::group_by(pw_id) %>%
    dplyr::arrange(met_id) %>%
    dplyr::mutate(met_ids = stringr::str_c(met_id, collapse = "|")) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(met_id,
              pw_id,
              pw_idx = dplyr::group_indices(., met_ids))

  pw_id_replacement <- row_data_indexed %>%
    dplyr::group_by(pw_idx) %>%
    dplyr::arrange(pw_id) %>%
    dplyr::mutate(tmp_ID = dplyr::first(pw_id),
           all_IDs = stringr::str_c(pw_id %>% unique(), collapse = "|")) %>%
    dplyr::ungroup() %>%
    dplyr::select(met_id,
                  ID = tmp_ID,
                  all_IDs) %>%
    dplyr::distinct() %>%
    dplyr::group_by(met_id) %>%
    dplyr::summarise(ID = stringr::str_c(ID, collapse = ", ")) %>%
    dplyr::mutate(ID = stringr::str_split(ID, ", ")) %>%
    dplyr::right_join(row_data, by = "met_id") %>%
    .$ID

  pwdb_summary_replacement <-
    dplyr::inner_join(metadata(D)$pathways[[pw_id]],
               dplyr::select(row_data_indexed, -met_id) %>%
                 dplyr::distinct(),
               by = c("ID" = "pw_id")) %>%
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
  rowData(D)[[pw_id]] <- pw_id_replacement
  metadata(D)$pathways[[pw_id]] <- pwdb_summary_replacement

  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('removed redundant pathway annotations using the %s column', pw_id)
    )

  D
}
