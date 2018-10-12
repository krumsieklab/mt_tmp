

#### source all files except myself... ugly but works for now ----
l=list.files(path=codes.makepath("packages/metabotools"),pattern="*.R",full.names=T)
l=l[!grepl('mtExamples.R',l)]
lapply(l, function(x){source(x,echo=F,verbose=F)})

source(codes.makepath('R/data/parseMetabolonFile.R'))


#### ----
D = parseMetabolonFile(codes.makepath("packages/metabotools/sampledata.xlsx"),"OrigScale",sumexp = T) %>% 
  mtFilter(metMax=0.2) %>% mtFilter(sampleMax=0.1) %>% mtQuotNorm()

assay(D)[1,1]
# metadata(D)
metadata(D)$preprocess[[1]]
metadata(D)$preprocess[[2]]

