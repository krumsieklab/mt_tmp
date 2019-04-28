library(tidyverse)
source(codes.makepath("MT/mt_internal_helpers.R"))

#' Remove redundant pathway annotations.
#' 
#' Remove identical pathways from an existing SummarizedExperiment
#' data structure that has a column of pathway annotations. 
#'
#' @param D \code{SummarizedExperiment} input
#' @param met_ID_col Column containing metabolite IDs
#' @param pw_col Column containing pathways IDs
#'
#' @return redundant pathway annotation from SE pw_col column will be filtered.
#'
#' @examples
#' # first annotate metabolites using smp_db and then remove redundant pathways
#' ... %>%
#'   mt_anno_pathways_HMDB(in_col = "HMDb_ID", out_col = "smp_db", 
#'   pwdb_name = "SMP", db_dir = codes.makepath("snippets/packages/metabotools_external/hmdb")) %>% 
#'   mt_anno_pathways_remove_redundant(met_ID_col = "HMDb_ID", pw_col = "smp_db") %>%
#' ...
#' 
#' @author Parviz Gomari
#' 

mt_anno_pathways_remove_redundant <- function(
  D,                  # SummarizedExperiment input
  met_ID_col,         # Column containing metabolite IDs
  pw_col              # Column containing pathways IDs
) {
  
  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  if(!met_ID_col %in% names(rowData(D)))
    stop(glue::glue("met_ID_col is invalid. please enter an existing column name."))
  
  if(!pw_col %in% names(rowData(D)))
    stop(glue::glue("pw_col is invalid. please enter an existing column name."))
  
  row_data <- 
    D %>% 
    rowData() %>% 
    as.data.frame() %>% 
    dplyr::select(met_ID = !!met_ID_col, 
                  pw_ID = !!pw_col) 
  
  mti_logstatus(glue::glue("creating grouping indices for the pathway IDs in {pw_col} "))
  row_data_indexed <- 
    row_data %>% 
    filter(!is.na(met_ID),
           pw_ID != "NULL") %>% 
    unnest(pw_ID) %>% 
    group_by(pw_ID) %>% 
    arrange(met_ID) %>% 
    mutate(met_IDs = str_c(met_ID, collapse = "|")) %>% 
    ungroup() %>% 
    transmute(met_ID, 
              pw_ID,
              pw_idx = group_indices(., met_IDs))
  
  pw_col_replacement <- row_data_indexed %>% 
    group_by(pw_idx) %>% 
    arrange(pw_ID) %>% 
    mutate(tmp_ID = dplyr::first(pw_ID),
           all_IDs = str_c(pw_ID %>% unique(), collapse = "|")) %>% 
    ungroup() %>% 
    dplyr::select(met_ID,
                  ID = tmp_ID, 
                  all_IDs) %>% 
    distinct() %>% 
    group_by(met_ID) %>%
    nest(ID, .key = ID) %>% 
    mutate(ID = 
             ID %>%
             unlist(recursive = FALSE) %>% 
             as.list()) %>% 
    right_join(row_data, by = "met_ID") %>% 
    .$ID
  
  pwdb_summary_replacement <- 
    inner_join(metadata(D)$pathways[[pw_col]],
               dplyr::select(row_data_indexed, -met_ID) %>% 
                 distinct(),
               by = c("ID" = "pw_ID")) %>% 
    group_by(pw_idx) %>% 
    arrange(ID) %>% 
    mutate(tmp_ID = dplyr::first(ID),
           all_IDs = str_c(ID, collapse = "|"),
           tmp_name = dplyr::first(pathway_name),
           all_pathway_names = str_c(pathway_name, collapse = "|")) %>% 
    filter(ID == tmp_ID) %>% 
    ungroup() %>% 
    dplyr::select(-c(tmp_ID, tmp_name, pw_idx))
  
  # replace rowData and metadata of D
  rowData(D)[[pw_col]] <- pw_col_replacement
  metadata(D)$pathways[[pw_col]] <- pwdb_summary_replacement
  
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('removed redundant pathway annotations using the %s column', pw_col)
    )
  
  D
}


if (FALSE) {
  # Example -----------------------------------------------------------------
  mt_logging(console=T) 
  D_alone <- 
    mt_files_load_metabolon(codes.makepath("MT/sampledata.xlsx"), "OrigScale") %>% 
    mt_anno_pathways_HMDB(in_col = "HMDb_ID", out_col = "smp_db", 
                         pwdb_name = "SMP", db_dir = codes.makepath("snippets/packages/metabotools_external/hmdb")) %>% 
    mt_anno_pathways_remove_redundant(met_ID_col = "HMDb_ID", pw_col = "smp_db")
}



