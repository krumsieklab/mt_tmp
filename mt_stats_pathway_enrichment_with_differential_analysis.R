#' Pathway enrichment using different methods.
#' 
#' This pathway enrichment analysis script calculats differential metabolites either
#' by a classical t-test or the GAGE package and then performs pathway enrichment
#' by using Fisher's exact test or GAGE.
#' 
#' Implemented approaches:
#' 1. GAGE: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-161
#' 2. Fisher's exact test. Will NOT scale() data before. Data matrix can have NAs.
#' 
#' @param D  \code{SummarizedExperiment} input
#' @param pw_col column containing pathways IDs
#' @param grouping_var column name of column containing grouping variables
#' @param control_grp_name name(s) of control group(s) found in grouping_var column
#' @param case_grp_name name(s) of case/phenotype group(s) found in grouping_var column
#' @param method one of: "gage", "fishers"
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
library(gage)


mt_stats_pathway_enrichment_with_differential_analysis <- function(
  D,          # SummarizedExperiment input
  pw_col,       # string, column name for pathway annotations
  grouping_var,
  control_grp_name,
  case_grp_name,
  method = "gage"
) {
  
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!method %in% c("fishers","gage")) stop("'method' must be either 'gage' or 'fishers'")
  
  meta_D <- metadata(D)
  
  if(!"pathways" %in% names(meta_D)) stop("'pathways' does not exist in current SummarizedExperiment input")
  
  # Check if given pathway column actually exists
  if (!pw_col %in% names(meta_D$pathways)) stop(sprintf("'%s' not found in metabolite annotations.", pw))
  
  pw_id_map <- 
    meta_D$pathways[[pw_col]] %>% 
    distinct(ID, pathway_name)
  
  row.data <- rowData(D) %>% as.data.frame()
  D.df <- assay(D) %>% as.data.frame()
  if(!("COMP_IDstr" %in% names(row.data))) {
    # If there isnt compound id (eg. if data from WCM core), make our own
    row.data$COMP_IDstr <- paste("cid", seq(nrow(row.data)), sep="")
    
    D.df$name <- rownames(D.df)
    nr <- nrow(D.df)
    
    
    row.data$new.name <- rownames(row.data)
    D.df <- inner_join(D.df, row.data %>% dplyr::select(new.name, COMP_IDstr),
                       by=c("name"="new.name"))
    row.data <- row.data %>% dplyr::select(-new.name)
    if(nr != nrow(D.df)) {
      stop("Number of rows changed after merge, something is wrong")
    }
    rownames(D.df) <- D.df$COMP_IDstr
    D.df <- D.df %>% dplyr::select(-c(COMP_IDstr, name))
  }
  
  
  geneset <- 
    row.data %>%
    dplyr::select(COMP_IDstr, !!rlang::sym(pw_col)) %>% 
    filter(!!rlang::sym(pw_col) != "NULL") %>% 
    unnest(!!rlang::sym(pw_col)) %>% 
    filter(!!rlang::sym(pw_col) != "NA") %>% 
    distinct()
  
  pw_names <- 
    geneset[[pw_col]] %>% 
    unique()
  
  geneset_list <- 
    sapply(
      pw_names, 
      simplify = FALSE,
      USE.NAMES = TRUE,
      function(pw_name) {
        geneset %>% 
          filter(!!rlang::sym(pw_col) == pw_name) %>% 
          .$COMP_IDstr
      }
    )
  
  # create indices for control and case groups
  # NOTE: it is assumed here that the order of sample groups in colData
  # follows the order of samples in the assay dataframe
  sample_indices <- 
    colData(D) %>% 
    as.data.frame() %>%
    dplyr::transmute(!!rlang::sym(grouping_var),
                     sample_idx = row_number())
  
  ctrl_samples <- 
    sample_indices %>% 
    filter(!!rlang::sym(grouping_var) %in% control_grp_name) %>% 
    .$sample_idx
  
  case_samples <- 
    sample_indices %>% 
    filter(!!rlang::sym(grouping_var) %in% case_grp_name) %>% 
    .$sample_idx
  
  
  if (method == "gage") {
    
    # perform pathway enrichment using the GAGE package
    enrichment_results <- 
      D.df %>%
      gage(gsets = geneset_list,
           ref = ctrl_samples,
           samp = case_samples,
           same.dir = FALSE, # suggsted for KEGG pathways
           compare = "unpaired", # unpaired samples
           set.size = c(5, 500)) %>% 
      .$greater %>% 
      as.data.frame() %>%
      rownames_to_column("pathway_ID") %>% 
      dplyr::select(pathway_ID,
                    p_value = p.val, 
                    p_value_adjusted = q.val,
                    mean_foldchange = stat.mean) %>% 
      filter(!is.na(p_value)) %>% 
      
      left_join(pw_id_map, by = c("pathway_ID" = "ID")) %>% 
      dplyr::select(pathway_name, everything()) 
    
    
  } else if (method == "fishers") {
    
    # Algorithm summary:
    # - calculate metabolite p-values per group
    # - adjust p-values using FDR
    # - assign significance to adjusted p-value at 0.05 level
    # - perform Fisher's exact test
    # - adjust Fisher's exact test p-value
    # - calculate mean fold change based on mean log value of cases and ctrls
    
    met_stats <- 
      D.df %>%
      t() %>% 
      as.data.frame() %>% 
      
      # assign sample index
      # NOTE: the assumption here is that each row is ordered 
      # as the rows in colData
      dplyr::mutate(sample_idx = row_number()) %>% 
      
      # select samples
      filter(sample_idx %in% c(ctrl_samples, case_samples)) %>% 
      dplyr::mutate(grp = if_else(sample_idx %in% ctrl_samples, "ctrl", "case")) %>% 
      dplyr::select(-sample_idx) %>% 
      
      # Calculate p-values and the group means
      summarise_if(is.numeric, list(p_value = ~ t.test(. ~ .data$grp)$p.value,
                                    mean_ctrl = ~ t.test(. ~ .data$grp)$estimate %>% as.numeric() %>% .[1],
                                    mean_case = ~ t.test(. ~ .data$grp)$estimate %>% as.numeric() %>% .[2])) %>% 
      gather(met_id, vals) %>% 
      separate(met_id, c("met_id", "val_type"), "_", extra = "merge") %>% 
      spread(val_type, vals) %>% 
      dplyr::mutate(p_value_adjusted = p.adjust(p_value, method = "fdr"))
    
    enrichment_results <- 
      met_stats %>% 
      
      # assign significance
      dplyr::mutate(significant = if_else(p_value_adjusted < 0.05, TRUE, FALSE),
                    n_total = n(),
                    n_total_sig = sum(significant)) %>% 
      inner_join(geneset, by = c("met_id" = "COMP_IDstr")) %>% 
      group_by(!!rlang::sym(pw_col)) %>% 
      
      # calculate summary numbers for Fisher's test
      dplyr::summarise(n_total = unique(n_total),
                       n_total_sig = unique(n_total_sig),
                       n_pw = n(),
                       n_pw_sig = sum(significant),
                       mean_case = mean(mean_case),
                       mean_ctrl = mean(mean_ctrl)) %>% 
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
                       mean_foldchange = mean_case - mean_ctrl) %>% 
      arrange(p_value)
    
  } else {
    
    stop("Bug.")
    
  }
  
  
  metadata(D)$pathways$enrichment_results <- 
    as_tibble(enrichment_results)
  
  
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('performed pathway enrichment on %s pathways using "%s"', 
                       nrow(enrichment_results), method)
    )
  
  D
}



if (FALSE) {
  
  # Example -----------------------------------------------------------------
  mt_logging(console=T) 
  D_pre <- 
    mt_files_load_metabolon(codes.makepath("MT/sampledata.xlsx"), "OrigScale") %>%
    mt_anno_pathways_HMDB(in_col = "HMDb_ID", out_col = "kegg_db", 
                          pwdb_name = "KEGG", 
                          #db_dir = codes.makepath("sharedcodes/packages/metabotools_external/hmdb/")) %>% 
                          db_dir = codes.makepath("snippets/packages/metabotools_external/hmdb/")) %>% 
    mt_anno_pathways_remove_redundant(met_ID_col = "HMDb_ID", pw_col = "kegg_db") %>% 
    mt_pre_filtermiss(metMax=0.2) %>%
    mt_pre_filtermiss(sampleMax=0.1) %>%
    # batch correction by variable BATCH_MOCK
    mt_pre_batch_median(batches = "BATCH_MOCK") %>%
    # quotient normalization
    mt_pre_norm_quot() %>%
    # logging
    mt_pre_trans_log() %>%
    # KNN imputation
    mt_pre_impute_knn()
  
  
  D_pw <- 
    D_pre %>% 
    mt_stats_pathway_enrichment_with_differential_analysis(
      "kegg_db", 
      grouping_var = "Group", 
      control_grp_name = "Vehicle", 
      case_grp_name = c("treatment2", "treatment1"),
      method = "gage"
    )
  
  metadata(D_pw)$pathways$enrichment_results
}

