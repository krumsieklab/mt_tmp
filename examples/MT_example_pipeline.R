#
# ------------ MetaboTools Example Pipeline ------------

# MetaboTools is an R statistical toolbox for performing Metabolomics analyses
# Built upon the packages SummarizedExperiment
# Designed to utilize the magrittr pipe operator for a smooth, continuous workflow
# Below is an example pipeline demonstrating the use of many of the functions provided by MetaboTools
# The MetaboTools user documentation provides an overview of all of the functions utilized in this example
# This example is divided into sections. Use SHIFT+CTRL+O to see the document outline.

library(MetaboTools)
library(tidyverse)
library(pathview)

zap()

# NTS: need to be able to load this via data()
file_data <- system.file("extdata", "example_data/simulated_data.xlsx", package = "MetaboTools")

# NOTE: IN THIS EXAMPLE WE DO NOT GROUP REUSABLE FUNCTION STEPS INTO FUNCTIONS BUT REPEAT THEM AS A DEMONSTRATION OF THE
# BROAD FUNCITON OF THE TOOLBOX

# PART 1 - STARTING A METABOTOOLS PIPELINE ----------------------------------------------------

D <-
  # validate checksum
  mt_files_checksum(file=file_data, checksum = "80afcd72481c6cf3dcf83342e3513699") %>%
  # load data - this function loads the assay data only
  #   alternative loading functions: mt_files_load_metabolon(), mt_files_load_metabolon_lipidomics(), mt_files_load_olink(),
  #     mt_files_load_UCD, mt_files_load_WCM, mt_files_load_nightingale, mt_files_load_metabolon_new_format()
  mt_files_data_xls(file=file_data, sheet="data", samples_in_row=T, ID_col="sample") %>%
  # load metabolite (rowData) annotations
  mt_files_anno_xls(file=file_data, sheet="metinfo",anno_type="metabolites", anno_ID="name", data_ID="name") %>%
  # load clinical (colData) annotations
  mt_files_anno_xls(file=file_data, sheet="clin", anno_type="samples", anno_ID="sample", data_ID="sample") %>%
  # # log assay dimensions and number of columns for both metabolite and clincial annotations
  mt_logging_datasetinfo() %>%
  # start timing
  mt_logging_tic() %>%
  {.}
# additional functions used at beginning of pipelines:
#   - mt_settings - set global settings for MetaboTools pipeline
#   - mt_flag_logged - for flagging a loaded dataset as already log transformed


# PART 2 - DATA CLEANING ----------------------------------------------------

D <- D %>%
  # heading
  mt_reporting_heading(strtitle = "Data Clean-up", lvl = 1) %>%
  # section text
  mt_reporting_text(text = "Filter samples that are missing values for Diagnosis,add sample annotation column with log10 of
                    PreBioPSA, convert sample annotaiton column Diagnosis to factors,filter metabolites that are missing values
                    for SUB_PATHWAY, log dataset information for this point of the pipeline.") %>%
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
#   - mt_pre_zeroToNA - for platforms that represent sub-LOD/missing values as zeros


# PART 3.1 - PREPROCESSING: FILTERING MISSING VALUES ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Preprocessing", lvl=1) %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Filtering", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot percent missingness for each metabolite before filtering, filter out metabolites with >= 50%
                    missingness, plot percent missingness for each metabolite after filtering, add missingness annotation
                    columns to both metabolite and sample annotation data frames.") %>%
  # plot missingness distribution
  mt_plots_qc_missingness(met_max=0.5) %>%
  # filter metabolites with more than 50% missing values per group
  mt_pre_filtermiss(met_max = 0.5, met_group = "Diagnosis") %>%
  # plot missingness distribution after filtering
  mt_plots_qc_missingness(met_max=0.5) %>%
  # add missingness percentage as annotation to samples (remaining missing)
  mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "metabolites", out_col = "missing") %>%
  {.}


# PART 3.2 - PREPROCESSING: NORMALIZATION ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Normalization", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot sample boxplots before normalization, apply median batch correction, perform quotient
                    normalization, plot boxplot with dilution factors from quotient normalization, plot sample boxplot after
                    normalization, log transform the data, impute missing data using knn, plot sample boxplot after imputation,
                    detect outliers, log dataset info, write pre-processed data to file.") %>%
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
  #   alternative imputation functions: mt_pre_impute_min
  mt_pre_impute_knn() %>%
  # plot sample boxplot after imputation
  mt_plots_sampleboxplot(color=Diagnosis, plottitle = "After imputation", logged = T) %>%
  # outlier detection (univariate)
  #   related function: mt_pre_outliercorrection
  # KELSEY: talk to Richa & Annalise
  mt_pre_outlier(method="univariate") %>%
  # correct metabolite abundances for Age
  #   alternative function: mt_pre_confounding_correction_stepwise_aic
  #mt_pre_confounding_correction(formula = ~ Age) %>%
  # print infos about dataset
  mt_logging_datasetinfo() %>%
  # write preprocessed data to Excel file
  #   other writing functions: mt_files_write_SE (write SummarizedExerpiment object)
  mt_files_write_xls(file = "PreprocessedData.xlsx") %>%
  {.}

# Additional pre-processing functions
#   - mt_pre_confounding_correction() - function for correcting confounding variables
#   - mt_pre_confounding_correction_stepwise_aic() - an alterenative function for correcting confounders that uses stepwise aic
# NOTES ON BEST PRACTICES: If incorporated in this pipeline, these functions would correct for the variable age such that
#     none of the following functions have to take care of those confounders anymore. If this function is included, confounders
#     should not be included in any of the following functions. It is generally agreed that including the confounders in the
#     linear models themselves is preferable to pre-correction.

# PART 4 - GET PATHWAY ANNOTATIONS ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Get Pathway Annotations", lvl = 1) %>%
  # get KEGG ids from HMDB ids
  mt_anno_metabolites_HMDB2KEGG(in_col = "HMDb", out_col = "KEGG_ids") %>%
  # get pathway annotations
  #   alternative functions: mt_anno_pathways_from_file, mt_anno_pathways_graphite, mt_anno_pathways_Uniprot
  mt_anno_pathways_HMDB(in_col = "HMDb", out_col = "pathway", pwdb_name = "KEGG") %>%
  # remove redundant
  mt_anno_pathways_remove_redundant(met_ID = "KEGG_ids", pw_ID = "pathway") %>%
  # write pathway annotations
  mt_files_write_pathwayannos(file="ExamplePipeline_PathwayAnnotations.xlsx", pwfield = "pathway") %>%
  {.}


# PART 5 - GLOBAL STATISTICS ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Global Statistics", lvl = 1) %>%
  # plot PCA
  mt_plots_PCA(scaledata = T, title = "scaled PCA - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP
  mt_plots_UMAP(scaledata = T, title = "scaled UMAP - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap
  mt_plots_pheatmap(scaledata = T, annotation_col = c("Diagnosis"), annotation_row = c("SUPER_PATHWAY"),
                    clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3) %>%
                    {.}


# PART 6.1 - STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS, METHOD: MISSINGNESS ANALYSIS ---------------------------------------

#create another SE object for first analysis branch (missingness)
D1 <- D

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Missingness analysis", lvl = 1) %>%
  # section text
  mt_reporting_text(text = "Perform missingness analysis to determine if NAs significantly accumulate in one of the Diagnosis
                    groups. Adjust output of test using multiple testing correction.") %>%
  # compute Fisher's exact test
  mt_stats_univ_missingness(comp_col="Diagnosis", stat_name="missingness") %>%
  # create p-value qq plot
  mt_plots_pvalqq(stat_name = "missingness") %>%
  # apply multiple testing correction
  #   alternative function: mt_post_multTest_effdim
  mt_post_multTest(stat_name="missingness", method="BH") %>%
  {.}


# PART 6.2 - STATISTICAL ANALYSIS, OUTCOME: AGE, METHOD: LINEAR REGRESSION -----------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Statistical Analysis", lvl = 1) %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    var = "Age",
                    stat_name = "Age met")%>%
  # create p-value qq plot
  mt_plots_pvalqq(stat_name = "Age met") %>%
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
  mt_plots_boxplot_scatter(stat_name = "Age met",
                           x = Age,
                           plot_type = "scatter",
                           metab_filter = p.adj < 1E-10, # made very small because otherwise would take an eternity to generate all plots
                           metab_sort = p.value,
                           annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
                           {.}


# PART 6.3 - STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS, METHOD: LINEAR REGRESSION (t-test) ----------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  #   alternative functions: mt_stats_univ_wilcox, mt_stats_univ_matrixeqtl
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
  mt_plots_boxplot_scatter(stat_name          ="Diagnosis met",
                           x                  = Diagnosis,
                           fill               = Diagnosis,
                           plot_type          = "box",
                           metab_filter       = p.adj < 0.05,
                           metab_sort         = p.value,
                           annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                           rows               = 3,
                           cols               = 3) %>%
                           {.}


# PART 7.1 - STATISTICAL RESULTS PRESENTATION: STATS BARPLOT & PATHVIEW ----------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Statisitcal Results Presentation", lvl = 1) %>%
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
  # NOTE: THIS FUNCTION CAN TAKE SEVERAL MINUTES TO RUN
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


# PART 7.2 - STATISTICAL RESULT PRESENTATION: MULTIPLE STATISTICS HEATMAP --------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(color_cutoff = 0.05) %>%
  {.}


# PART 7.3 - STATISTICAL RESULT PRESENTATION: RESULT COMPARISON ----------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Result Comparison", lvl = 2) %>%
  # comparison plot
  mt_plots_compare2stats(stat1 = "Age met", filter1 = p.adj < 0.05,
                         D2 = D1, stat2 = "Diagnosis met", filter2 = p.adj < 0.05,
                         filterop = "OR") %>%
                         {.}




# PART 8 - PARTIAL CORRELATION NETWORK ----------------------------------------------------

# THE DATA TABLE IS TOO LARGE TO INCLUDE IN AN HTML FILE

#D1 <- D1 %>%
# heading for html file
#  mt_reporting_heading(strtitle = "Partial Correlation Network", lvl = 2) %>%
# compute partial correlation matrix
#  mt_stats_multiv_net_GeneNet(stat_name = "GGM") %>%
# add multiple testing correction
#  mt_post_multTest(stat_name = "GGM", method = "BH") %>%
# plot network and color according to age analysis
#  mt_plots_net(stat_name = "GGM", corr_filter = p.adj < 0.05, node_coloring = "Age met") %>%
#  {.}


# PART 9 - PATHWAY AGGREGATION ANALYSIS ----------------------------------------------------
# This is now first aggregating the metabolite matrix into pathways creating a new matrix of pathway
# conentration values, and then repeating the parts of the same pipeline as above
# In a real scenario, you would include a statistical analysis performed on metabolites AND have to pick between part 9 and the sections above

# create another SE object for third analysis branch (pathway)
D2 <- D

D2 <- D2 %>%
  # aggregate metabolites in the same pathways
  mt_modify_aggPW(pw_col = "pathway", method = "aggmean") %>%

  # STATISTICAL ANALYSIS, OUTCOME: AGE
  # heading for html file
  mt_reporting_heading(strtitle = "Pathway Aggregation Analysis", lvl = 1) %>%
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
  mt_plots_boxplot_scatter(stat_name = "Age pw",
                           x = Age,
                           plot_type = "scatter",
                           metab_filter = p.adj < 1E-10,
                           metab_sort = p.value,
                           annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}",
                           rows = 3,
                           cols = 3) %>%

  # STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS
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
  mt_plots_boxplot_scatter(stat_name          ="Diagnosis pw",
                           x                  = Diagnosis,
                           fill               = Diagnosis,
                           plot_type          = "box",
                           metab_filter       = p.adj < 0.05,
                           metab_sort         = p.value,
                           annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}")

D2 <- D2 %>%
  # MULTIPLE STATISTICS HEATMAP
  # heading for html file
  mt_reporting_heading(strtitle = "Statistical Results Presentation", lvl = 2) %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(color_cutoff = 0.05) %>%

  # COMPARE STATISTICAL RESULTS
  # heading for html file
  mt_reporting_heading(strtitle = "Result Comparison", lvl = 2) %>%
  # comparison plot from pw analysis
  mt_plots_compare2stats(stat1 = "Age pw", filter1 = p.adj < 0.05,
                         D2 = D2, stat2 = "Diagnosis pw", filter2 = p.adj < 0.05,
                         filterop = "OR")


# PART 10 - SUB PATHWAY ANALYSIS ----------------------------------------------------
# create another SE object for third analysis branch (pathway)
D3 <- D

D3 <- D3 %>%
  # aggregate metabolites in the same pathways
  mt_modify_aggPW(pw_col = "SUB_PATHWAY", method = "aggmean") %>%

  # STATISTICAL ANALYSIS, OUTCOME: AGE
  # heading for html file
  mt_reporting_heading(strtitle = "Sub Pathway Aggregation Analysis", lvl = 1) %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    var = "Age",
                    stat_name = "Age sub pw") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "Age sub pw", method = "BH") %>%
  # add stats logging
  mt_logging_statsinfo(stat_name = "Age sub pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pvalhist(stat_names = "Age sub pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Age sub pw",
                   x = statistic,
                   metab_filter = p.adj < 0.05,
                   colour = p.adj < 0.05) %>%
  # scatter plot
  mt_plots_boxplot_scatter(stat_name = "Age sub pw",
                           x = Age,
                           plot_type = "scatter",
                           metab_filter = p.adj < 1E-10,
                           metab_sort = p.value,
                           annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}",
                           rows = 3,
                           cols = 3) %>%

  # STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS
  # heading for html file
  mt_reporting_heading(strtitle = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Diagnosis,
                   stat_name = "Diagnosis sub pw") %>%
  # add fold change
  mt_post_addFC(stat_name = "Diagnosis sub pw") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "Diagnosis sub pw", method = "BH") %>%
  # add stats logging
  mt_logging_statsinfo(stat_name = "Diagnosis sub pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pvalhist(stat_names = "Diagnosis sub pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Diagnosis sub pw",
                   x = fc,
                   metab_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_boxplot_scatter(stat_name          ="Diagnosis sub pw",
                           x                  = Diagnosis,
                           fill               = Diagnosis,
                           plot_type          = "box",
                           metab_filter       = p.adj < 0.05,
                           metab_sort         = p.value,
                           annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}")

D3 <- D3 %>%
  # MULTIPLE STATISTICS HEATMAP
  # heading for html file
  mt_reporting_heading(strtitle = "Statistical Results Presentation", lvl = 2) %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(color_cutoff = 0.05) %>%

  # COMPARE STATISTICAL RESULTS
  # heading for html file
  mt_reporting_heading(strtitle = "Result Comparison", lvl = 2) %>%
  # comparison plot from pw analysis
  mt_plots_compare2stats(stat1 = "Age sub pw", filter1 = p.adj < 0.05,
                         D2 = D3, stat2 = "Diagnosis sub pw", filter2 = p.adj < 0.05,
                         filterop = "OR") %>%
  {.}




# PART 11 - SUB PATHWAY ANALYSIS  ----------------------------------------------------
# create another SE object for third analysis branch (pathway)
D4 <- D

D4 <- D4 %>%
  # aggregate metabolites in the same pathways
  mt_modify_aggPW(pw_col = "SUPER_PATHWAY", method = "aggmean") %>%

  # STATISTICAL ANALYSIS, OUTCOME: AGE
  # heading for html file
  mt_reporting_heading(strtitle = "Super Pathway Aggregation Analysis", lvl = 1) %>%
  # add tag
  mt_reporting_tag(tag_name = "Beginning of Super PW Age Analysis") %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Age analysis", lvl = 2) %>%
  # Pearson correlation
  mt_stats_univ_cor(method = "pearson",
                    var = "Age",
                    stat_name = "Age super pw") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "Age super pw", method = "BH") %>%
  # add stats logging
  mt_logging_statsinfo(stat_name = "Age super pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pvalhist(stat_names = "Age super pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Age super pw",
                   x = statistic,
                   metab_filter = p.adj < 0.05,
                   colour = p.adj < 0.05) %>%
  # scatter plot
  mt_plots_boxplot_scatter(stat_name = "Age super pw",
                           x = Age,
                           plot_type = "scatter",
                           metab_filter = p.adj < 1E-10,
                           metab_sort = p.value,
                           annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}",
                           rows = 3,
                           cols = 3) %>%

  # STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS
  # add tag
  mt_reporting_tag(tag_name = "Beginning of Super PW Diagnosis Analysis") %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Diagnosis analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Diagnosis,
                   stat_name = "Diagnosis super pw") %>%
  # add fold change
  mt_post_addFC(stat_name = "Diagnosis super pw") %>%
  # add multiple testing correction
  mt_post_multTest(stat_name = "Diagnosis super pw", method = "BH") %>%
  # add stats logging
  mt_logging_statsinfo(stat_name = "Diagnosis super pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pvalhist(stat_names = "Diagnosis super pw") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "Diagnosis super pw",
                   x = fc,
                   metab_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_boxplot_scatter(stat_name          ="Diagnosis super pw",
                           x                  = Diagnosis,
                           fill               = Diagnosis,
                           plot_type          = "box",
                           metab_filter       = p.adj < 0.05,
                           metab_sort         = p.value,
                           annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}")

D4 <- D4 %>%
  # MULTIPLE STATISTICS HEATMAP
  # heading for html file
  mt_reporting_heading(strtitle = "Statistical Results Presentation", lvl = 2) %>%
  # heading for html file
  mt_reporting_heading(strtitle = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(color_cutoff = 0.05) %>%

  # COMPARE STATISTICAL RESULTS
  # heading for html file
  mt_reporting_heading(strtitle = "Result Comparison", lvl = 2) %>%
  # comparison plot from pw analysis
  mt_plots_compare2stats(stat1 = "Age super pw", filter1 = p.adj < 0.05,
                         D2 = D4, stat2 = "Diagnosis super pw", filter2 = p.adj < 0.05,
                         filterop = "OR") %>%
  {.}

# PART 12 - COMPARE SUPER- AND SUB- PATHWAY ANALYSES  ----------------------------------------------------
D4 <- D4 %>% mt_plots_equalizer(comp1 = "Age super pw",
                                D2 = D3,
                                comp2 = "Age sub pw",
                                legend.fine="sub pathway",
                                legend.coarse='super pathway',
                                vertline.fine = p.adj < 0.1,
                                vertline.coarse = p.adj < 0.1) %>%
  # end timing
  mt_logging_toc() %>%
  {.}

# PART 13 - CREATE ANALYSIS REPORTS ----------------------------------------------------

# OPTIONAL - this function will remove all plot entries in the pipeline
#   - mt_stripresults(strip = "plots")

# metabolite analysis html report
D1 %>% mt_reporting_html(outfile = "Example_Pipeline_Metabolite_Analysis.html",
                         title = "Example Pipeline - Statistical Analysis")
# pathway analysis html report
D2 %>% mt_reporting_html(outfile = "Example_Pipeline_Pathway_Analysis.html",
                         title = "Example Pipeline - Pathway Aggregation Analysis")

# sub-pathway analysis html report
D3 %>% mt_reporting_html(outfile = "Example_Pipeline_Sub_Pathway_Analysis.html",
                         title = "Example Pipeline - Sub Pathway Aggregation Analysis")

# super-pathway analysis html report
D4 %>% mt_reporting_html(outfile = "Example_Pipeline_Super_Pathway_Analysis.html",
                         title = "Example Pipeline - Super Pathway Aggregation Analysis",
                         start.after = "Beginning of Super PW Age Analysis")

# create a combined report with both analsyses
mt_reporting_html_nonLinear(pipelines = list(D1, D2, D3, D4), outfile = "ExamplePipeline_CombinedReport.hmtl",
                            title = "Combined Report")
