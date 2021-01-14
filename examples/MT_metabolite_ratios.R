### -- Stand-Alone Example: Replace Metabolites with Ratios -- ###
# This script demonstrates the functionality of mt_modify_ratios() and mt_post_pgain().

library(MetaboTools)

D <-
  # load data
  mt_files_load_metabolon(codes.makepath("Mt/sampledata.xlsx"), "OrigScale") %>%
  # timing start
  mt_logging_tic() %>%

  ### Preprocessing ---------
  # heading
  mt_reporting_heading("Preprocessing") %>%
  mt_reporting_heading("Part 1", lvl=2) %>%
  # sample boxplot
  mt_plots_sampleboxplot() %>%
  # missingness plot
  mt_plots_qc_missingness() %>%
  # filter metabolites with >20% missing values, then samples with >10% missing values
  mt_pre_filtermiss(met_max=0.2) %>%
  mt_pre_filtermiss(sample_max=0.1) %>%
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batches = "BATCH_MOCK") %>%
  # heading
  mt_reporting_heading("Part 2", lvl=2) %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  # check if there is any correlation between normalization factors and outcomes (bad sign if so)
  mt_plots_qc_dilutionplot(comp="num1") %>%
  mt_plots_qc_dilutionplot(comp="Group") %>%
  # logging
  mt_pre_trans_log() %>%
  # KNN imputation
  mt_pre_impute_knn() %>%



  # Modify data frame to represent metabolites as ratios ---------
  # GGM
  mt_stats_multiv_net_GeneNet(stat_name ="GGM") %>%
  mt_post_multTest(stat_name = "GGM", method="fdr") %>%
  mt_modify_ratios(nbr_stat_name = "GGM", nbr_edge_filter = p.adj<0.5, nbr_neighborhood = 2) %>% # liberal cutoff because this is a mock dataset

  # Perform analysis on metabolite ratios ---------
  mt_reporting_heading("Statistics") %>%
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
  mt_post_multTest(stat_name = "comp", method = "BH") %>%
  # p-value histogram
  mt_plots_pvalhist(stat_name = "comp") %>%
  # pgains
  mt_post_pgain(stat_name = "comp") %>%
  # Volcano plot as overview of results
  mt_plots_volcano(stat_name     = "comp",
                   metab_filter = p.adj < 0.1 & pgain > 2,
                   colour       = p.adj < 0.1 & pgain > 2) %>%

                   {.}


# Generate HTML report ---------

D %>% mt_reporting_html(outfile="example_ratios.html", output.calls = T)

