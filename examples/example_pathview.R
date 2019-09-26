
#### load libraries ----

zap()
# load MT 
mt.quickload()

#### part 1, preprocess data ----

mt_logging(console=T) 
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
  mt_pre_impute_knn()

#### part 2, run differential analysis ----

D <- D %>%
  ###### </SAME AS example_simplepipeline.R>
  ### heading
  mt_reporting_heading("Statistics") %>%
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group, 
    samplefilter = (Group %in% c("treatment1","treatment2")),
    name         = "comp",
    mc.cores     = 1
  ) %>%
  # add fold change
  mt_post_addFC(statname = "comp") %>%
  # add multiple testing correction
  mt_post_multTest(statname = "comp", method = "BH") %>%
  # p-value histogram
  mt_plots_pvalhist() %>%
  # Volcano plot as overview of results
  mt_plots_volcano(statname     = "comp",
                   x = statistic,
                   metab_filter = p.adj < 0.1,
                   colour       = p.value < 0.05)

#### part 3, create pathview plots ----

D <- D %>%
  # add column to rowData with KEGG identifiers (computed from HMDB)
  mt_anno_metabolites_HMDB2KEGG(col_input = "HMDb_ID", col_output = "KEGG_identifiers")

D <- D %>%
  # plot all metabolites in the top 10 most frequent pathway annotations
  mt_plots_pathview(met.id="KEGG_identifiers",
                    n.pathways = 3,
                    # kegg pathway files will be created in a folder called "Pathview_database" inside the current working directory
                    path.database = "./Pathview_database",
                    # output will be created in a folder called "Pathview_output1" inside the current working directory
                    path.output = "./Pathview_output1")
  # plot all metabolites in the top 10 most frequent pathway annotations
  # use the results of the statistical analysis "comp" to color and filter metabolites
  
D <- D %>%
  mt_plots_pathview(met.id="KEGG_identifiers",
                    n.pathways = 3,
                    # take results from statistical analysis called "comp"
                    statname = "comp",
                    # color scale function
                    color.scale = -sign(fc)*log10(p.value),
                    # set color range
                    color.range = 3.5,
                    # metabolite filtering condition
                    metab.filter = p.value < 0.05,
                    # get pathway list only from filtered metabolites
                    show.only.filtered = TRUE,
                    # kegg pathway files will be created in a folder called "Pathview_database" inside the current working directory
                    path.database = "./Pathview_database",
                    # output will be created in a folder called "Pathview_output2" inside the current working directory
                    path.output = "./Pathview_output2",
                    # set to false to speed-up (output files will be bigger in size)
                    same.layer = FALSE)

D <- D %>%
  mt_plots_pathview(met.id="KEGG_identifiers",
                    n.pathways = 3,
                    # take results from statistical analysis called "comp"
                    statname = "comp",
                    # color scale function
                    color.scale = -sign(fc)*log10(p.value),
                    # set color range
                    color.range = 3.5,
                    # metabolite filtering condition
                    metab.filter = p.value < 0.05,
                    # get pathway list from all metabolites
                    show.only.filtered = FALSE,
                    # kegg pathway files will be created in a folder called "Pathview_database" inside the current working directory
                    path.database = "./Pathview_database",
                    # output will be created in a folder called "Pathview_output2" inside the current working directory
                    path.output = "./Pathview_output2",
                    # set to false to speed-up (output files will be bigger in size)
                    same.layer = FALSE)
