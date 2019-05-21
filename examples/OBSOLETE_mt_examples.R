####
####
#### OLD SANDBOX FOR EXAMPLE CODE
####
####


#### source all files except myself... ugly but works for now ----
require(tidyverse)
zap()  # to be safe, erases all vars, can be added to .Rprofile.....   zap <- function(){lst<-ls(envir=.GlobalEnv); lst<-lst[!(lst %in% c("zap","codes.makepath","store","restore","debugstore"))]; rm(list=lst, envir=.GlobalEnv) }

source(codes.makepath("MT/quickload.R"))

#### run single pipeline, preprocessing, linear models, multiple testing correction etc. no networks ----
mt_logging(console=T) 
D <- 
  # load data
  mt_files_load_metabolon(codes.makepath("Mt/sampledata.xlsx"), "OrigScale") %>%
  # timing start
  mt_logging_tic() %>% 
  # PCA, colored by some rowData() fields... this function shows 3 PCs
  mt_plots_PCA_mult(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>%
  # sample boxplot
  mt_plots_sampleboxplot() %>%
  # missingness plot
  mt_plots_qc_missingness() %>%
  # filter metabolites with >20% missing values, then samples with >10% missing values
  mt_pre_filtermiss(metMax=0.2) %>%
  mt_pre_filtermiss(sampleMax=0.1) %>%
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batches = "BATCH_MOCK") %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  # check if there is any correlation between normalization factors and outcomes (bad sign if so)
  mt_plots_qc_dilutionplot(comp="num1") %>%
  mt_plots_qc_dilutionplot(comp="Group") %>%
  # logging
  mt_pre_trans_log() %>%
  # KNN imputatio
  mt_pre_impute_knn() %>%
  # outlier detection (univariate)
  mt_pre_outlier(method="univariate", thr=4, perc=0.5) %>%
  # final sample boxplot
  mt_plots_sampleboxplot(color=Group, plottitle = 'final') %>%
  # PCA, colored by some rowData() fields... this function shows 2 PCs
  mt_plots_PCA(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>%
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group, 
    samplefilter = (Group %in% c("Li_2","Li_5")),
    name         = "Li's",
    mc.cores     = 1
  ) %>%
  # add fold changes to result tables
  mt_post_addFC(statname = "Li's") %>%
  # add multiple testing correction
  mt_post_multTest(statname = "Li's", method = "BH") %>%
  # p-value histogram
  mt_plots_pvalhist() %>%
  # Volcano plot as overview of results
  mt_plots_volcano(statname     = "Li's",
                   metab_filter = p.adj < 0.1,
                   colour       = p.value < 0.05) %>%
  # boxplots
  mt_plots_boxplot(statname           = "Li's",
                   x                  = Group,
                   fill               = Group,
                   correct_confounder = ~BRADFORD_PROTEIN + BATCH_MOCK,
                   metab_filter       = p.value<0.01,
                   metab_sort         = p.value,
                   annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
                   rows               = 2,
                   cols               = 2) %>%
  mt_stats_multiv_net_GeneNet(name = "pcor") %>%
  mt_post_multTest(statname ="pcor", method = "BH") %>%
  mt_plots_net(statname = "pcor", corr_filter = p.value<0.01, node_coloring = "Li's") %>%
  # final timing
  mt_logging_toc()

# show all functions that have been run
metadata(D)$results %>% map('fun') %>% map_chr(str_c, collapse = " - ")

# plot all plots to PDF file
pdf("output.pdf")
walk(D %>% mti_res_get_plots(), plot)
dev.off()
# alternative code for the same result
D %>% mti_plot_all_tofile(file="output.pdf")

# plot all plots to R plot window (will stack up all plots to browse through)
walk(D %>% mti_res_get_plots(), plot)


# # missing value analysis
# D_missing <-
#   mt_files_load_metabolon(codes.makepath("MT/sampledata.xlsx"), "OrigScale") %>%
#   mt_stats_univ_missingness(comp = 'Group', name='miss')
# 
# mti_get_stat_by_name(D_missing, "miss")

# pathways
# D_sub <- D_alone %>% 
#   mt_modify_aggPW(pw="SUB_PATHWAY", method="aggmean") %>% 
#   mt_stats_univ_lm(
#     formula      = ~ Group, 
#     samplefilter = (Group %in% c("Li_2","Li_5")),
#     name         = "Li's SUB"
#   ) %>%
#   mt_plots_equalizer(comp1="Li's SUB", D2=D_alone, comp2="Li's", legend.fine="metabolite", legend.coarse='sub pathway')
# 
# # walk(D_sub %>% mti_res_get_plots(), plot)
# metadata(D_sub)$results[[27]]
# 
# D_super <- D_sub %>% 
#   mt_modify_aggPW(pw="SUPER_PATHWAY", method="aggmean") %>% 
#   mt_stats_univ_lm(
#     formula      = ~ Group, 
#     samplefilter = (Group %in% c("Li_2","Li_5")),
#     name         = "Li's SUPER"
#   ) %>%
#   mt_plots_equalizer(comp1="Li's SUPER", D2=D_sub, comp2="Li's SUB", legend.fine="sub pathway", legend.coarse='super pathway')
# 
# # walk(D_super %>% mti_res_get_plots(), plot)
# metadata(D_super)$results[[30]]

# mt_stats_multiv_net_GeneNet(name="GeneNetpcor") %>%
#   mt_post_multTest(statname = "GeneNetpcor", method = "BH") %>%
#   mt_plots_net(statname = "GeneNetpcor", corr_filter = p.adj < 0.5, node_coloring="Li's") %>%  
#   
# #### run pathway analysis with HMDB KEGG annotations ----
# D_kegg <- D_alone %>% 
#   mt_anno_pathways_HMDB(in_col = "HMDb_ID", out_col = "kegg_db", 
#                        pwdb_name = "KEGG", db_dir = codes.makepath("packages/metabotools_external/hmdb")) %>% 
#   mt_modify_aggPW(pw="kegg_db", method="aggmean") %>% 
#   mt_stats_univ_lm(
#     formula      = ~ Group, 
#     samplefilter = (Group %in% c("Li_2","Li_5")),
#     name         = "Li's KEGG"
#   ) %>% 
#   mt_post_addFC(statname = "Li's KEGG") %>%
#   mt_post_multTest(statname = "Li's KEGG", method = "BH") %>%
#   mt_plots_volcano(statname     = "Li's KEGG",
#                    metab_filter = p.adj < 0.2,
#                    colour       = p.value < 0.05)
# 
# # last plot
# r <- metadata(D_kegg)$results
# r[[length(r)]]$output
# 
# 
# 
# #### integrate all results into structure ----
# 
# r <- MTResultCollector$new()
# 
# # r$add(D_alone)
# # r$add(D_sub)
# r$addMultiple(D_alone, D_sub, D_super)
# 
# # play around with graph walking
# root <- r$graph_roots()
# r$graph_next(root)
# # branching point from alone to pathways
# branch <- metadata(D_alone)$results[[length(metadata(D_alone)$results)]]$uuid
# r$graph_next(branch)
# 
# igraph.options(plot.layout=layout.auto, vertex.size=5, label.cex=0.1 )
# # igraph.options(plot.layout=layout_with_sugiyama, vertex.size=5)
# 
# plot(simplify(r$graph))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### run preprocessing ----
# message("\nPREPROCESS")
# mt_logging(console=T) 
# D <- 
#   mt_files_load_metabolon(codes.makepath("MT/sampledata.xlsx"), "OrigScale") %>%
#   mt_plots_PCA_mult(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>%
#   mt_plots_sampleboxplot() %>%
#   mt_plots_qc_missingness() %>%
#   mt_pre_filtermiss(metMax=0.2) %>%
#   mt_pre_filtermiss(sampleMax=0.1) %>%
#   mt_pre_batch_median(batches = "BATCH_MOCK") %>%
#   mt_pre_norm_quot() %>%
#   mt_pre_trans_log() %>%
#   mt_plots_qc_dilutionplot(comp="num1") %>%
#   mt_plots_qc_dilutionplot(comp="Group") %>%
#   mt_pre_impute_knn() %>%
#   mt_plots_sampleboxplot(color=Group) %>%
#   mt_plots_PCA(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) 
# 
# 
# 
# 
# 
# #### create analysis pipeline (note %>>>% operator)
# message("\nCREATE PIPELINE")
# require(gestalt)
# my_sexy_analysis <- {.} %>>>%
#   mt_stats_univ_lm(
#     formula      = ~ Group, 
#     samplefilter = (Group %in% c("Li_2","Li_5")),
#     name         = "Li's",
#     mc.cores     = 1
#   ) %>>>%
#   mt_stats_univ_lm(
#     formula      = ~ num2 + Group, 
#     name         = "num",
#     mc.cores     = 1
#   ) %>>>%
#   mt_post_multTest(statname = "Li's", method = "BH") %>>>%
#   mt_post_addFC(statname = "Li's") %>>>%
#   mt_plots_volcano(statname     = "Li's",
#                    metab_filter = p.adj < 0.1,
#                    colour       = p.value < 0.05) %>>>%
#   mt_plots_volcano(statname     = "Li's",
#                    x            = -rank(abs(statistic)),
#                    metab_filter = rank(p.value) <= 10,
#                    colour       = SUPER_PATHWAY) %>>>%
#   mt_internal_void(param=0.5) %>>>%
#   mt_plots_boxplot(statname           = "Li's",
#                    x                  = Group,
#                    fill               = Group,
#                    correct_confounder = ~BRADFORD_PROTEIN + BATCH_MOCK,
#                    metab_filter       = 1:n() < 6,  
#                    metab_sort         = p.value,
#                    annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
#                    rows               = 2,
#                    cols               = 2) %>>>%
#   mt_plots_pvalhist(.) %>>>%  ## THE DOT IS TECHNICALLY REDUNDANT BUT NECESSARY SO THAT FUNCTION NAME IS PRESERVED
#   mt_internal_void(param=2)
# 
# 
# #### apply pipeline to individual metabolites and ratios ----
# message("\nSINGLE")
# D_single <- D %>%
#   my_sexy_analysis() 
# metadata(D_single)$results %>% map('fun') %>% map_chr(str_c, collapse = " - ")
# 
# D_ratio <- D %>%
#   mt_modify_ratios() %>%
#   mt_modify_filter_metabolites(metab_filter = 1:n() < 1000) %>%
#   mt_modify_filter_samples(sample_filter = Group %in% c("Li_2","Li_5")) %>%
#   my_sexy_analysis() %>%
#   mt_post_pgain(statname = "Li's")
# 
# #### print results table (incl. p-gain and mult test correction)
# message("\nRATIOS")
# D_ratio %>% mti_get_stat_by_name("Li's") %>% arrange(desc(pgain), p.adj) %>% head()
# 
# 
# #### apply pipeline to pathways ----
# message("\nPATHWAYS")
# D_sub <- D %>%
#   mt_modify_aggPW(pw="SUB_PATHWAY", method="aggmean") %>%
#   my_sexy_analysis() 
# ## rowData(D_sub)
# 
# 
# 
# 
# #### plotting ----
# 
# #### extract plots of all analyses ----
# message("\nEXTRACT PLOTS")
# pl_qc     <- D %>% mti_res_get_plots(unlist=T)
# pl_single <- D_single %>% mti_res_get_plots(unlist=T) %>% tail(-length(pl_qc))
# pl_ratio  <- D_ratio  %>% mti_res_get_plots(unlist=T) %>% tail(-length(pl_qc))
# pl_pw  <- D_sub  %>% mti_res_get_plots(unlist=T) %>% tail(-length(pl_qc))
# 
# #### all plots of all analyses ----
# pdf("output_all.pdf")
# walk(pl_qc    , plot)
# walk(pl_single, plot)
# walk(pl_ratio , plot)
# walk(pl_pw , plot)
# dev.off()
# 
# #### just all plots of one analysis ----
# pdf("output_single.pdf")
# walk( D_single %>% mti_res_get_plots(unlist=T)    , plot)
# dev.off()
# 
