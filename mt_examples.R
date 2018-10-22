

#### source all files except myself... ugly but works for now ----
zap()  # to be safe, erases all vars, can be added to .Rprofile.....   zap <- function(){lst<-ls(envir=.GlobalEnv); lst<-lst[!(lst %in% c("zap","sysdiab.makepath","store","restore","debugstore"))]; rm(list=lst, envir=.GlobalEnv) }

l=list.files(path=codes.makepath("packages/metabotools"),pattern="*.R$",full.names=T)
l=l[!grepl('*examples*',l)]
lapply(l, function(x){source(x,echo=F,verbose=F)})




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
  mt_stats_univ_glm(formula = num1 ~ M + num2)




#### plotting ----

# get all plots
pl <- D %>% mti_res_get_plots()

# plot in R window
sapply(pl,plot)

# export
pdf("output.pdf")
sapply(pl, plot)
dev.off()



#### show statistics ----




