### -- MT Average Duplicate Samples & Metabolites -- ###
# MetaboTools provides functions for calculating CV when QC samples are provided

library(MetaboTools)
data("D_ex1")

# D_ex1 has two types of QC pools
colData(D_ex1)$QCID %>% table()

# Add cv from both QC pools to rowData
D_ex1 <- D_ex1 %>%
  mt_modify_cv(QCID == "QC_1", col_lab = "QC1_cv") %>%
  mt_modify_cv(QCID == "QC_2", col_lab = "QC2_cv")
