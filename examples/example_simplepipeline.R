#####
##### MetaboTools pipeline
#####

#### run pipeline ----

# load MT
mt.quickload()

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
  mt_pre_filtermiss(metMax=0.2) %>%
  mt_pre_filtermiss(sampleMax=0.1) %>%
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
  # outlier detection (multivariate) & visualization
  # mt_pre_outlier(method="mahalanobis", pval=0.01, reduce.dim = T) %>%
  mt_plots_PCA(color='outlier_mahalanobis') %>%
  # final sample boxplot
  mt_plots_sampleboxplot(color=Group, plottitle = 'final') %>%
  # PCA, colored by some rowData() fields... this function shows 2 PCs
  mt_plots_PCA(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>%
  
  ###
  mt_reporting_heading("Statistics") %>%
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group, 
    samplefilter = (Group %in% c("treatment1","treatment2")),
    name         = "comp",
    mc.cores     = 1
  ) %>%
  # add fold changes to result tables
  mt_post_addFC(statname = "comp") %>%
  # add multiple testing correction
  mt_post_multTest(statname = "comp", method = "BH") %>%
  # p-value histogram
  mt_plots_pvalhist() %>%
  # Volcano plot as overview of results
  mt_plots_volcano(statname     = "comp",
                   metab_filter = p.adj < 0.1,
                   colour       = p.value < 0.05) %>%
  
  ###
  # heading
  mt_reporting_heading("All boxplots") %>%
  # boxplots
  mt_plots_boxplot(statname           = "comp",
                   x                  = Group,
                   fill               = Group,
                   correct_confounder = ~BATCH_MOCK,
                   metab_filter       = p.value<0.01,
                   metab_sort         = p.value,
                   annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
                   rows               = 2,
                   cols               = 2) %>%
  # final timing
  mt_logging_toc() %>%
  
  # testing void (should not occur)
  mt_internal_void()




#### generate and knit markdown ----

D %>% mt_reporting_quickhtml(outfile="example_simplepipeline.html", output.calls = T)

