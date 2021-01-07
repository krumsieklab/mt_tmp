### -- MT Pathway Enrichment Analysis -- ###
# Perform pathway enrichment analysis using the function mt_stats_pathway_enrichment()


library(MetaboTools)
# without differential analysis

# Preprocess dataset ------------------------------------------------------
D_pre <-
  mt_files_load_metabolon(file=codes.makepath("Mt/sampledata.xlsx"), sheet="OrigScale") %>%
  mt_anno_pathways_HMDB(in_col = "HMDb_ID",
                        out_col = "kegg_db",
                        pwdb_name = "KEGG",
                        db_dir = data.makepath("MT_precalc/hmdb/")) %>%
  mt_anno_pathways_remove_redundant(met_ID = "HMDb_ID", pw_ID = "kegg_db") %>%
  mt_pre_filtermiss(met_max=0.2) %>%
  mt_pre_filtermiss(sample_max=0.1) %>%
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batches = "BATCH_MOCK") %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  # logging
  mt_pre_trans_log() %>%
  # KNN imputation
  mt_pre_impute_knn() %>%
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group,
    sample_filter = (Group %in% c("treatment1","treatment2")),
    stat_name         = "comp",
    n_cores     = 1
  ) %>%
  # add fold changes to result tables
  mt_post_addFC(stat_name = "comp") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "comp", method = "BH")




# Apply mt_stats_pathway_enrichment ---------------------------------------

D_pw <-
  D_pre %>%
  mt_stats_pathway_enrichment(pw_col = "kegg_db",
                              stat_name = "comp",
                              cutoff = 0.4)


# show results
metadata(D_pw)$pathways$enrichment_results
