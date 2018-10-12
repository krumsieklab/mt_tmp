
#### ----
D = parseMetabolonFile(codes.makepath("packages/metabotools/sampledata.xlsx"),"OrigScale",sumexp = T) %>% 
  mtFilter(metMax=0.2) %>% mtFilter(sampleMax=0.1) %>% mtQuotNorm()

assay(D)[1,1]
# metadata(D)
metadata(D)$preprocess[[1]]
metadata(D)$preprocess[[2]]

