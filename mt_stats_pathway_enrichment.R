#' Pathway enrichment script using statistical analysis results from metabotools pipeline.
#' 
#' This script runs a pathway enrichment analysis using Fisher's exact test from statistics
#' calculated by mt_stats_univ_lm().
#' 
#' Implemented approaches:
#' 1. Fisher's exact test. 
#' 
#' @param D  \code{SummarizedExperiment} input
#' @param pw_col column containing pathways IDs
#' @param stat_name name of statistical analysis to use with mti_get_stat_by_name()
#' @param cutoff cutoff to use for assigning whether a metabolite/gene is significant. Used in Fisher's exact test
#'
#' @return $pathways$enrichment_results: a dataframe containing the pathway enrichment results
#' 
#' @examples
#' %>% mt_stats_pathway_enrichment("kegg_db", 
#'                                 grouping_var = "Group", 
#'                                 control_grp_name = "Vehicle", 
#'                                 case_grp_name = c("treatment2", "treatment1") %>%
#' 
#' 
#' @author Parviz Gomari
#' 

library(tidyverse)


mt_stats_pathway_enrichment <- function(
  D,          # SummarizedExperiment input
  pw_col,      # string, column name for pathway annotations
  stat_name,
  cutoff = 0.05
) {
  
  stopifnot("SummarizedExperiment" %in% class(D))
  
  meta_D <- metadata(D)
  
  if(!"pathways" %in% names(meta_D)) stop("'pathways' does not exist in current SummarizedExperiment input")
  
  # Check if given pathway column actually exists
  if (!pw_col %in% names(meta_D$pathways)) stop(sprintf("'%s' not found in metabolite annotations.", pw))
  
  # have a check for wheter stat_name exists in D?
  
  pw_id_map <- 
    meta_D$pathways[[pw_col]] %>% 
    distinct(ID, pathway_name)
  
  
  geneset <- 
    rowData(D) %>% 
    as.data.frame() %>% 
    dplyr::select(COMP_IDstr, !!rlang::sym(pw_col)) %>% 
    filter(!!rlang::sym(pw_col) != "NULL") %>% 
    unnest(!!rlang::sym(pw_col)) %>% 
    filter(!!rlang::sym(pw_col) != "NA") %>% 
    distinct()
  
  
  # Algorithm summary:
  # - calculate metabolite p-values per group
  # - adjust p-values using FDR
  # - assign significance to adjusted p-value at 0.05 level
  # - perform Fisher's exact test
  # - adjust Fisher's exact test p-value
  # - calculate mean fold change based on mean log value of cases and ctrls
  
  enrichment_results <- 
    mti_get_stat_by_name(D, stat_name) %>% 
    
    # assign significance
    dplyr::mutate(significant = if_else(p.adj < cutoff, TRUE, FALSE),
                  n_total = n(),
                  n_total_sig = sum(significant)) %>% 
    inner_join(geneset, by = c("var" = "COMP_IDstr")) %>% 
    group_by(!!rlang::sym(pw_col)) %>% 
    
    # calculate summary numbers for Fisher's test
    dplyr::summarise(n_total = unique(n_total),
                     n_total_sig = unique(n_total_sig),
                     n_pw = n(),
                     n_pw_sig = sum(significant),
                     mean_fc = mean(fc)) %>% 
    filter(n_pw >= 5) %>% 
    
    # calculate contingency table entries
    dplyr::mutate(s_p = n_pw_sig,
                  ns_p = n_pw - n_pw_sig,
                  s_np = n_total_sig - n_pw_sig,
                  ns_np = n_total - (s_p + ns_p + s_np)) %>% 
    rowwise() %>% 
    dplyr::mutate(p_value = 
                    fisher.test(matrix(c(s_p, s_np, ns_p, ns_np), nrow = 2)) %>% 
                    .$p.value) %>% 
    ungroup() %>% 
    dplyr::rename(ID = !!rlang::sym(pw_col)) %>% 
    left_join(pw_id_map, by = "ID") %>% 
    dplyr::transmute(pathway_name, 
                     pathway_ID = ID,
                     p_value, 
                     p_value_adjusted = p.adjust(p_value, method = "fdr"), 
                     mean_foldchange = mean_fc) %>% 
    arrange(p_value)
  
  metadata(D)$pathways$enrichment_results <- 
    as_tibble(enrichment_results)
  
  
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("performed pathway enrichment on %s pathways using Fihser's exact test", 
                       nrow(enrichment_results))
    )
  
  D
}