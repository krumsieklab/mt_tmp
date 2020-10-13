#####
##### MetaboTools Example Pipeline
#####
# MetaboTools is an R statistical toolbox for performing Metabolomics analyses
# Built upon the packages SummarizedExperiment
# Designed to utilize the magrittr pipe operator for a smooth, continuous workflow
# Below is an example pipeline demonstrating the use of many of the functions provided by MetaboTools
# The MetaboTools user documentation provides an extensive overview of all of the functions utilized in this example

library(MetaboTools)
library(tidyverse)

zap()

# NTS: DELETE THIS WHEN EXAMPLE FINISHED
setwd("/Users/kelsey/Desktop/KrumsiekLab/")

# NTS: need to be able to load this via data()
file_data <- "simulated_data.xlsx"

##############################################################################################################################
# PART 1 - STARTING A METABOTOOLS PIPELINE
##############################################################################################################################
D <-
  # load data - this function loads the assay data only
  #   alternative loading functions: mt_files_load_metabolon, mt_files_load_metabolon_lipidomics, mt_files_load_olink,
  #     mt_files_load_UCD, mt_files_load_WCM
  mt_files_data_xls(file=file_data, sheet="data", samples_in_row=T, ID_col="sample") %>%
  # validate checksum
  mt_files_checksum(file=file_data, checksum = "80afcd72481c6cf3dcf83342e3513699") %>%
  # load metabolite (rowData) annotations
  mt_files_anno_xls(file=file_data, sheet="metinfo",anno_type="metabolites", anno_ID="name", data_ID="name") %>%
  # load clinical (colData) annotations
  mt_files_anno_xls(file=file_data, sheet="clin", anno_type="samples", anno_ID="sample", data_ID="sample") %>%
  # # log assay dimensions and number of columns for both metabolite and clincial annotations
  mt_logging_datasetinfo() %>%
  # start timing
  mt_logging_tic() %>%
  {.}
  # additional functions commonly used at start of pipeline:
  #   - mt_flag_logged - allows user to flag if data is already logged

##############################################################################################################################
# PART 2 - DATA CLEANING
##############################################################################################################################

D <- D %>%
  # heading
  mt_reporting_heading(strtitle = "Data Clean-up", lvl = 1) %>%
  # filter samples
  mt_modify_filter_samples(sample_filter = !is.na(Diagnosis)) %>%
  # create additional variable
  mt_modify_mutate(anno_type = "samples", col_name = "PreBioPSALog", term = log10(PreBioPSA)) %>%
  # modify variable to factor
  mt_modify_applytoanno(anno_type = "samples", col_name = "Diagnosis", fun = as.factor) %>%
  # remove metabolites with no pathway annotation
  mt_modify_filter_metabolites(metab_filter = !is.na(SUB_PATHWAY)) %>%
  # log assay dimensions and number of columns for both metabolite and clinical annotations
  mt_logging_datasetinfo() %>%
  {.}
  # additional data cleaning function:
  #   - mt_pre_zeroToNA - for platforms that represent sub-LOB/missing values as zeros



##############################################################################################################################
# PART 3.1 - PREPROCESSING: MISSING ANALYSIS
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Preprocessing", lvl = 1) %>%

  # heading for html file
  mt_reporting_heading(strtitle = "Missingness Analysis", lvl = 2) %>%
  # compute Fisher's exact test
  mt_stats_univ_missingness(comp_col="Diagnosis", stat_name="missingness") %>%
  # apply multiple testing correction
  #   alternative function: mt_post_multTest_effdim
  mt_post_multTest(stat_name="missingness", method="BH") %>%
  {.}

##############################################################################################################################
# PART 3.2 - PREPROCESSING: FILTERING
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Filtering", lvl = 2) %>%
  # plot missingness distribution
  mt_plots_qc_missingness(met_max=0.5) %>%
  # filter metabolites with more than 50% missing values per group
  mt_pre_filtermiss(met_max = 0.5, met_group = "Diagnosis") %>%
  # plot missingness distribution after filtering
  mt_plots_qc_missingness(met_max=0.5) %>%
  # add missingness percentage as annotation to samples
  mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "metabolites", out_col = "missing") %>%
  {.}

##############################################################################################################################
# PART 3.3 - PREPROCESSING: NORMALIZATION
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Normalization", lvl = 2) %>%
  # plot sample boxplots
  mt_plots_sampleboxplot(color=Diagnosis, plottitle = "Original", logged = T) %>%
  # apply batch correction
  #   alternative batch correction function: mt_pre_batch_COMBAT
  mt_pre_batch_median(batches = "BOX.NUMBER") %>%
  # plot sample boxplots after batch correction
  mt_plots_sampleboxplot(color=Diagnosis, plottitle = "After batch correction", logged = T) %>%
  # normalize abundances using probabilistic quotient
  #   alternative normalization function: mt_pre_norm_external
  mt_pre_norm_quot(met_max = 0.2, ref_samples = Diagnosis==0) %>%
  # show dilution plot
  mt_plots_qc_dilutionplot(comp="Diagnosis") %>%
  # plot sample boxplots after normalization
  mt_plots_sampleboxplot(color=Diagnosis, plottitle = "After normalization", logged = T) %>%
  # log transform
  #   other data transformation functions: mt_pre_trans_exp, mt_pre_trans_relative, mt_pre_trans_scale
  mt_pre_trans_log() %>%
  # impute missing values using knn
  #   alternative imputation function: mt_pre_impute_knn_multicore
  mt_pre_impute_knn() %>%
  # plot sample boxplot after imputation
  mt_plots_sampleboxplot(color=Diagnosis, plottitle = "After imputation", logged = T) %>%
  # outlier detection (univariate)
  #   related function: mt_pre_outliercorrection
  mt_pre_outlier(method="univariate") %>%
  # print infos about dataset
  mt_logging_datasetinfo() %>%
  # write preprocessed data to Excel file
  #   other writing functions: mt_files_write_SE (write SummarizedExerpiment object)
  mt_files_write_xls(file = "PreprocessedData.xlsx") %>%
  {.}

##############################################################################################################################
# PART 3.3 - PREPROCESSING: GLOBAL STATISTICS
##############################################################################################################################

# NTS: potentially start a branch here - leaving these plots out could save some space

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Global Statistics", lvl = 2) %>%
  # plot PCA
  mt_plots_PCA(scaledata = T, title = "scaled PCA - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP
  mt_plots_UMAP(scaledata = T, title = "scaled UMAP - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap
  mt_plots_pheatmap(scaledata = T, annotation_col = c("Diagnosis"), annotation_row = c("SUPER_PATHWAY"),
                    clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3) %>%
  {.}

##############################################################################################################################
# PART 4.1 - METABOLITE ANALYSIS: AGE ANALYSIS
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Metabolite Analysis", lvl = 1) %>%

  # heading for html file
  mt_reporting_heading(strtitle = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    var = "Age",
                    stat_name = "Age met")%>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "Age met", method = "BH") %>%
  # add stats logging
  mt_logging_statsinfo(stat_name = "Age met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_files_write_stats(file = "AgeAnalysis.xlsx", metname = "BOCHEMICAL") %>%
  # pvalue histogram
  mt_plots_pvalhist(stat_names = "Age met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Age met",
                   x = statistic,
                   metab_filter = p.adj < 0.05,
                   colour = p.adj < 0.05) %>%
  # scatter plot
  mt_plots_scatter(stat_name = "Age met",
                   x = Age,
                   metab_filter = p.adj < 1E-10, # made very small because otherwise would take an eternity to generate all plots
                   metab_sort = p.value,
                   annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}",
                   rows = 3,
                   cols = 3) %>%
  # correct metabolite abundances for Age
  #   alternative function: mt_pre_confounding_correction_stepwise_aic
  mt_pre_confounding_correction(formula = ~ Age) %>%
  {.}

##############################################################################################################################
# PART 4.2 - METABOLITE ANALYSIS: DIAGNOSIS ANALYSIS
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  #   alternative function: mt_stats_univ_wilcox
  mt_stats_univ_lm(formula = ~ Diagnosis,
                   stat_name = "Diagnosis met") %>%
  # add fold change
  mt_post_addFC(stat_name = "Diagnosis met") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "Diagnosis met", method = "BH") %>%
  # add stats logging
  mt_logging_statsinfo(stat_name = "Diagnosis met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_files_write_stats(file = "DiagnosisAnalysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pvalhist(stat_names = "Diagnosis met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Diagnosis met",
                   x = fc,
                   metab_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_boxplot(stat_name           ="Diagnosis met",
                   x                  = Diagnosis,
                   fill               = Diagnosis,
                   metab_filter       = p.adj < 0.05,
                   metab_sort         = p.value,
                   annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                   rows               = 3,
                   cols               = 3) %>%
  {.}

##############################################################################################################################
# PART 4.3 - METABOLITE ANALYSIS: BAR PLOT
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Barplot", lvl = 2) %>%
  # create statsbarplots
  mt_plots_statsbarplot(stat_name = c("Age met", "Diagnosis met"),
                        metab_filter = p.adj < 0.05,
                        aggregate = "SUB_PATHWAY",
                        colorby = "SUPER_PATHWAY",
                        yscale = "count",
                        sort = T,
                        assoc_sign = "statistic") %>%
  mt_plots_statsbarplot(stat_name = c("Age met", "Diagnosis met"),
                        metab_filter = p.adj < 0.05,
                        aggregate = "SUB_PATHWAY",
                        colorby = "SUPER_PATHWAY",
                        yscale = "fraction",
                        sort = T,
                        assoc_sign = "statistic") %>%
  {.}

##############################################################################################################################
# PART 4.4 - METABOLITE ANALYSIS: PARTIAL CORRELATION NETWORK
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Partial Correlation Network", lvl = 2) %>%
  # compute partial correlation matrix
  mt_stats_multiv_net_GeneNet(stat_name = "GGM") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "GGM", method = "BH") %>%
  # plot network and color according to age analysis
  mt_plots_net(stat_name = "GGM", corr_filter = p.adj < 0.05, node_coloring = "Age met") %>%
  {.}

##############################################################################################################################
# PART 5.1 - PATHWAY ANALYSIS: GET ANNOTATIONS
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Pathway Analysis", lvl = 1) %>%

  # heading for html file
  mt_reporting_heading(strtitle = "Get Annotations", lvl = 2) %>%
  # get KEGG ids from HMDB ids
  mt_anno_metabolites_HMDB2KEGG(in_col = "HMDb", out_col = "KEGG_ids") %>%
  # get pathway annotations
  #   alternative functions: mt_anno_pathways_from_file, mt_anno_pathways_graphite, mt_anno_pathways_Uniprot
  mt_anno_pathways_HMDB(in_col = "HMDb", out_col = "pathway", pwdb_name = "KEGG") %>%
  # remove redundant
  mt_anno_pathways_remove_redundant(met_ID = "KEGG_ids", pw_ID = "pathway") %>%
  # write pathway annotations
  mt_files_write_pathwayannos()
  {.}

##############################################################################################################################
# PART 5.2 - PATHWAY ANALYSIS: PATHVIEW PLOT
##############################################################################################################################

D <- D %>%
  mt_plots_pathview(met.id="KEGG_ids",
                    species = "hsa",
                    # n.pathways = 5,
                    # take results from statistical analysis called "Age met"
                    stat_name = "Age met",
                    # color scale function
                    color.scale = sign(statistic),
                    # metabolite filtering condition
                    metab.filter = p.adj < 0.05,
                    # get pathway list only from filtered metabolites
                    show.only.filtered = TRUE,
                    # kegg pathway files will be created in a folder called "Pathview_database" inside the current working directory
                    path.database = "./Pathview_database",
                    # output will be created in a folder called "Pathview_output" inside the current working directory
                    path.output = "./Pathview_output",
                    # set to false to speed-up (output files will be bigger in size)
                    same.layer = FALSE,
                    add.pwname.suffix = TRUE) %>%
  {.}

##############################################################################################################################
# PART 5.3 - PATHWAY ANALYSIS: AGE ANALYSIS
##############################################################################################################################

D <- D %>%
  # aggregate metabolites in the same pathways
  mt_modify_aggPW(pw_col = "pathway", method = "aggmean") %>%

  # heading for html file
  mt_reporting_heading(strtitle = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    var = "Age",
                    stat_name = "Age pw") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "Age pw", method = "BH") %>%
  # add stats logging
  mt_logging_statsinfo(stat_name = "Age pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pvalhist(stat_names = "Age pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Age pw",
                   x = statistic,
                   metab_filter = p.adj < 0.05,
                   colour = p.adj < 0.05) %>%
  # scatter plot
  mt_plots_scatter(stat_name = "Age pw",
                   x = Age,
                   metab_filter = p.adj < 1E-10,
                   metab_sort = p.value,
                   annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}",
                   rows = 3,
                   cols = 3) %>%
  {.}

##############################################################################################################################
# PART 5.4 - PATHWAY ANALYSIS: DIAGNOSIS ANALYSIS
##############################################################################################################################

# heading for html file
mt_reporting_heading(strtitle = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Diagnosis,
                   stat_name = "Diagnosis pw") %>%
  # add fold change
  mt_post_addFC(stat_name = "Diagnosis pw") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "Diagnosis pw", method = "BH") %>%
  # add stats logging
  mt_logging_statsinfo(stat_name = "Diagnosis pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pvalhist(stat_names = "Diagnosis pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Diagnosis pw",
                   x = fc,
                   metab_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_boxplot(stat_name          ="Diagnosis pw",
                   x                  = Diagnosis,
                   fill               = Diagnosis,
                   metab_filter       = p.adj < 0.05,
                   metab_sort         = p.value,
                   annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                   rows               = 3,
                   cols               = 3) %>%
  {.}

##############################################################################################################################
# PART 6.1 - RESULT COMPARISON: METABOLITE
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Result Comparison", lvl = 1) %>%

  # heading for html file
  mt_reporting_heading(strtitle = "Metabolite", lvl = 2) %>%
  # comparison plot
  mt_plots_compare2stats(stat1 = "Age met", filter1 = p.adj < 0.05,
                         D2 = D, stat2 = "Diagnosis met", filter2 = p.adj < 0.05,
                         filterop = "OR") %>%
  {.}

##############################################################################################################################
# PART 6.2 - RESULT COMPARISON: PATHWAY
##############################################################################################################################

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Pathway", lvl = 2) %>%

  # comparison plot from pw analysis
  mt_plots_compare2stats(stat1 = "Age pw", filter1 = p.adj < 0.05,
                         D2 = D, stat2 = "Diagnosis pw", filter2 = p.adj < 0.05,
                         filterop = "OR") %>%
  # end timing
  mt_logging_toc()
  {.}

##############################################################################################################################
# PART 7 - CREATE HTML DOCUMENT
##############################################################################################################################

# create html
D %>% mt_reporting_html(outfile = "ExamplePipeline.html", title = "Example Pipeline")
