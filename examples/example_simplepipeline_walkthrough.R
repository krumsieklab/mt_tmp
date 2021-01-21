#####
##### MetaboTools pipeline
#####

#### run pipeline ----

# load MT
library(MetaboTools)

D <-
  # load data
  mt_load_files_metabolon(file = codes.makepath("Mt/sampledata.xlsx"), sheet = "OrigScale") %>%
  ## new SE
  ## assay: 244 rows and 12 columns; column names are 1:12, row names are from metabolon "BIOCHEMICAL" column
  ## rowData: 244 rows with 16 metabolite annotation columns
  ## colData:12 rows with 9 sample annotation columns
  ## metadata:
  ##  result: standard info; no output or output2
  ##  sessionInfo: session information at time of run
  ##  parseInfo: file name and sheet name

  # timing start
  mt_logging_tic() %>%
  ## metadata:
  ##  result: standard info

  ###
  # heading
  mt_reporting_heading("Preprocessing") %>%
  ## metadata:
  ##  result: standard info; output - list with level and title string
  mt_reporting_heading("Part 1", lvl=2) %>%
  ## metadata:
  ##  result: standard info; output - list with level and title string
  # sample boxplot
  mt_plots_sampleboxplot() %>%
  ## metadata:
  ##  result: standard info; output - list of plots "p"
  # missingness plot
  mt_plots_qc_missingness() %>%
  ## metadata
  ##  result: standard info; output - list of plots "plots"
  # filter metabolites with >20% missing values, then samples with >10% missing values
  mt_pre_filtermiss(met_max=0.2) %>%
  ## assay: filter out rows with missingness >= 20% (244 -> 186)
  ## rowData: filter out rows removed from assay (244 -> 186)
  ## metadata:
  ##  result: standard info; output - list of rows kept
  mt_pre_filtermiss(sample_max=0.1) %>%
  ## assay: filter out columns with missingness >= 10% (12 -> 12)
  ## colData: filter out rows removed from assay (12 -> 12)
  ## metadata:
  ##  result: standard info; output - list of columns kept
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batches = "BATCH_MOCK") %>%
  ##
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
  # heatmap
  mt_plots_pheatmap(scaledata = T) %>%

  ###
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
  mt_plots_pvalhist() %>%
  # Volcano plot as overview of results
  mt_plots_volcano(stat_name     = "comp",
                   metab_filter = p.adj < 0.1,
                   colour       = p.value < 0.05) %>%

  ###
  # heading
  mt_reporting_heading("All boxplots") %>%
  # boxplots
  mt_plots_boxplot(stat_name           = "comp",
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
  {.}




#### generate and knit markdown ----

D %>% mt_reporting_html(outfile="/Users/kelsey/Desktop/example_simplepipeline_walkthrough.html", output.calls = T)

