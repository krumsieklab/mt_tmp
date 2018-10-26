
#### source all files except myself... ugly but works for now ----
zap()  # to be safe, erases all vars, can be added to .Rprofile.....   zap <- function(){lst<-ls(envir=.GlobalEnv); lst<-lst[!(lst %in% c("zap","codes.makepath","store","restore","debugstore"))]; rm(list=lst, envir=.GlobalEnv) }

l=list.files(path=codes.makepath("packages/metabotools"),pattern="*.R$",full.names=T)
l=l[!grepl('*examples*',l)]
lapply(l, function(x){source(x,echo=F,verbose=F)})

# library(operators)


#### run ----
mt_logging(console=T) 
D <- 
  mt_files_load_metabolon(codes.makepath("packages/metabotools/sampledata.xlsx"), "OrigScale") %>%
  mt_plots_PCA_mult(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>%
  mt_plots_sampleboxplot() %>%
  mt_plots_qc_missingness() %>%
  mt_pre_filtermiss(metMax=0.2) %>%
  mt_pre_filtermiss(sampleMax=0.1) %>%
  mt_pre_trans_log() %>%
  mt_pre_norm_quot() %>%
  mt_plots_qc_dilutionplot(comp="num1") %>%
  mt_plots_qc_dilutionplot(comp="Group") %>%
  mt_pre_impute_knn() %>%
  mt_plots_sampleboxplot(color=Group) %>%
  mt_plots_PCA(color=Group, shape=BATCH_MOCK, size=NUM_MOCK)


#### create analysis pipeline (note %>>>% operator)
require(gestalt)
my_sexy_analysis <- {.} %>>>%
    mt_stats_univ_lm(
        formula      = ~ Group, 
        samplefilter = (Group %in% c("Li_2","Li_5")),
        name         = "Li's",
        mc.cores     = 2
    ) %>>>%
    mt_stats_univ_lm(
        formula      = ~ num2 + Group, 
        name         = "num",
        mc.cores     = 2
    ) %>>>%
    mt_post_multTest(statname = "Li's", method = "BH") %>>>%
    mt_post_addFC(statname = "Li's") %>>>%
    mt_plots_volcano(statname     = "Li's",
                     metab_filter = p.adj < 0.1,
                     colour       = p.value < 0.05) %>>>%
    mt_plots_volcano(statname     = "Li's",
                     x            = -rank(abs(statistic)),
                     metab_filter = rank(p.value) <= 10,
                     colour       = SUPER_PATHWAY) %>>>%
    mt_plots_boxplot(statname           = "Li's",
                     x                  = Group,
                     fill               = Group,
                     correct_confounder = ~BRADFORD_PROTEIN + BATCH_MOCK,
                     metab_filter       = 1:n() < 6,  
                     metab_sort         = p.value,
                     annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
                     rows               = 2,
                     cols               = 2)

#### apply pipeline to individual metabolites and ratios
D_single <- D %>%
    my_sexy_analysis() 

D_ratio <- D %>%
    mt_modify_ratios() %>%
    mt_modify_filter_metabolites(metab_filter = 1:n() < 1000) %>%
    mt_modify_filter_samples(sample_filter = Group %in% c("Li_2","Li_5")) %>%
    my_sexy_analysis() %>%
    mt_post_pgain(statname = "Li's")

#### print results table (incl. p-gain and mult test correction)
D_ratio %>% mti_get_stat_by_name("Li's") %>% arrange(desc(pgain), p.adj)

#### plotting ----

  
# get all plots
pl_qc     <- D %>% mti_res_get_plots(unlist=T)
pl_single <- D_single %>% mti_res_get_plots(unlist=T) %>% tail(-length(pl_qc))
pl_ratio  <- D_ratio  %>% mti_res_get_plots(unlist=T) %>% tail(-length(pl_qc))

# plot in R window
sapply(pl_qc,plot)

# export
pdf("output.pdf")
sapply(pl_qc    , plot)
sapply(pl_single, plot)
sapply(pl_ratio , plot)
dev.off()


#### show statistics ----
