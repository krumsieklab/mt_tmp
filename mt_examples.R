

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
  mt_plots_sampleboxplot() %>%
  mt_pre_filtermiss(metMax=0.2) %>%
  mt_pre_filtermiss(sampleMax=0.1) %>%
  mt_pre_trans_log() %>%
  mt_pre_norm_quot() %>%
  mt_pre_impute_knn() %>%
  mt_plots_sampleboxplot(color=Group) %>%
  mt_plots_PCA(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>% 
  mt_stats_univ_lm(
    formula      = ~ Group, 
    samplefilter = (Group %in% c("Li_2","Li_5")),
    name         = "Li's"
    ) %>%
    mt_plots_boxplot(statname           = "Li's",
                   x                  = Group,
                   fill               = Group,
                   correct_confounder = ~BRADFORD_PROTEIN + BATCH_MOCK,
                   metab_filter       = p.value < 0.05 & !is.na(SUPER_PATHWAY),  
                   metab_sort         = p.value,
                   annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
                   rows               = 2,
                   cols               = 2)




#### plotting ----

# get all plots
pl <- D %>% mti_res_get_plots(unlist=T)

# plot in R window
sapply(pl,plot)

# export
pdf("output.pdf")
sapply(pl, plot)
dev.off()


#### show statistics ----




