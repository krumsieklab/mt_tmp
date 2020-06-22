#####
##### MetaboTools pipeline
#####

#### network-based ratios (NBR) example

#### run pipeline ----

# load MT
library(MetaboTools)

D <-
  # load data
  mt_files_load_metabolon(codes.makepath("Mt/sampledata.xlsx"), "OrigScale") %>%
  # timing start
  mt_logging_tic() %>%
  
  ###
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
  
  
  
  # ratios
  
  # GGM
  mt_stats_multiv_net_GeneNet(stat_name ="GGM") %>%
  mt_post_multTest(stat_name = "GGM", method="fdr") %>%
  
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group,
    sample_filter = (Group %in% c("treatment1","treatment2")),
    stat_name         = "comp",
    n_cores     = 1
  )  %>%
  
  # plot network
  mt_plots_net(stat_name = "GGM", node_coloring = "comp")
  
  {.}

# D%>% MetaboTools:::mti_get_stat_by_name("comp") %>% View()


#### generate and knit markdown ----

D %>% mt_reporting_html(outfile="example_net.html", output.calls = T)

