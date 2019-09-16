
#### load libraries ----

library(AnnotationDbi)
library(org.Hs.eg.db)
library(graphite)

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
  # linear model on multiple groups, i.e. ANOVA
  mt_stats_univ_lm(
    formula      = ~ Group, 
    name         = "comp",
    mc.cores     = 1
  ) %>%
  # add multiple testing correction
  mt_post_multTest(statname = "comp", method = "BH") %>%
  # p-value histogram
  mt_plots_pvalhist() %>%
  # Volcano plot as overview of results
  mt_plots_volcano(statname     = "comp",
                   x = statistic,
                   metab_filter = p.adj < 0.1,
                   colour       = p.value < 0.05)


#### part 3, gnerate KEGG pathway annotations ----

# download kegg pathways
pwdb <- pathways(species = "hsapiens", database = "kegg")
pwdb <- mcmapply(function(pwname) {
  pw <- pwdb[[pwname]]
  pw %>% 
    graphite::edges(which = "mixed") %>% # "metabolites", "proteins", or "mixed"
    mutate(name = pwname,
           ID = pathwayId(pw))
},
names(pwdb),
SIMPLIFY = F,
mc.cleanup = T,
mc.cores = 3)

# build just one big dataframe with all pathway informations
pwdf <- do.call(rbind, pwdb)

# find metabolite pathway annotations
m_anno <- lapply(rowData(D)$KEGG[!is.na(rowData(D)$KEGG)], function(x) {
  pwdf$ID[pwdf$dest==x] %>% unique()
  
})
names(m_anno) <- rowData(D)$KEGG[!is.na(rowData(D)$KEGG)]

# build one long list
m_anno_list <- do.call(c, m_anno)

# find most common pathway for metabolites
pw_met <- m_anno_list %>% table() %>% as.data.frame() 
colnames(pw_met) <- c("pathway","Freq")
# pathway list ordered according to the number of metabolites with that annotation
pw_met <- pw_met[order(pw_met$Freq,decreasing = TRUE),]
# remove ":" from pathway ids for pathview
pw_met$pathway <- gsub(":", "", pw_met$pathway)

#### part 4, create pathview plots ----

D <- D %>%
  # plot all metabolites and all genes in the top 5 most frequent pathway annotations for metabolites
  mt_plots_pathview(met.id="KEGG", 
                    pathway.id = pw_met$pathway[1:5], 
                    # plots will be created in a folder called "Pathview_met" inside the current working directory
                    kegg.dir = "./Pathview_met") %>%
  # plot only significant metabolites in the top 5 most frequent pathway annotations for metabolites
  mt_plots_pathview(met.id="KEGG", 
                    pathway.id = pw_met$pathway[1:5], 
                    statname = "comp", 
                    # only metabolites with p.adj < 0.05 in the statistitical comparison called "comp" will be included in the plots
                    metab_filter = p.adj < 0.05,
                    # plots will be created in a folder called "Pathview_met_sing" inside the current working directory
                    kegg.dir = "./Pathview_met_sign")
