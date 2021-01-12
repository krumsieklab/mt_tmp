### -- MT Average Duplicate Samples & Metabolites -- ###
# MetaboTools provides functions for averaging duplicate samples or metabolites

library(MetaboTools)

# Average Dupcliate Samples -----------------
# Load toy SE object with duplicated samples
data("D_ex1")

# 100 samples and 100 metabolites
D_ds %>% assay() %>% as.data.frame() %>% dim()
# 5 samples are duplicates according to Sample_ID
colData(D_ds)$Sample_ID %>% duplicated() %>% sum()

# use MetaboTools to average + combine samples
D_ds <- D_ds %>%
  mt_modify_averagesample(group_by = "Sample_ID")

# 95 unique samples and 100 metabolites
D_ds %>% assay() %>% as.data.frame() %>% dim()
colData(D_ds)$Sample_ID %>% duplicated() %>% sum()


# Average Duplciate Metabolites -----------------
# Load toy SE object with duplicated metabolites
data("D_ex2")

D_ex2 %>% assay() %>% as.data.frame() %>% dim()
rowData(D_ex2)$COMP_ID %>% duplicated() %>% sum()

D_ex2 <- D_ex2 %>%
  mt_modify_averagemetabolite(avg_by = "COMP_ID")

D_ex2 %>% assay() %>% as.data.frame() %>% dim()
rowData(D_ex2)$COMP_ID %>% duplicated() %>% sum()

