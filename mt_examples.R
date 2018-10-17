

#### source all files except myself... ugly but works for now ----
l=list.files(path=codes.makepath("packages/metabotools"),pattern="*.R",full.names=T)
l=l[!grepl('*examples*',l)]
lapply(l, function(x){source(x,echo=F,verbose=F)})

source(codes.makepath('R/data/parseMetabolonFile.R'))


#### run ----
D = parseMetabolonFile(codes.makepath("packages/metabotools/sampledata.xlsx"),"OrigScale",sumexp = T) %>% 
  mt_plots_sampleboxplot(colorby="Group") %>% 
  mt_pre_filtermiss(metMax=0.2) %>%
  mt_pre_filtermiss(sampleMax=0.1) %>%
  mt_pre_trans_log() %>% 
  mt_pre_norm_quot() %>% 
  mt_pre_impute_knn() %>% 
  mt_plots_sampleboxplot(colorby="Group") %>%
  mt_plots_PCA(colorby="Group",shapeby="BATCH_MOCK")



#### ----

# show last plot
metadata(D)$plots[[length(metadata(D)$plots)]]


# plot all to PDF
pdf("output.pdf")
sapply(metadata(D)$plots, plot)
dev.off()

