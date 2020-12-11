### mt_modify_averagesample Example ###

load("MT_sample_dataset.rda")

D <- toy_dataset %>%
  # Add cv from both QC pools to colData
  mt_modify_cv(QCID == "QC_1", col_lab = "QC1_cv")%>%
  mt_modify_cv(QCID == "QC_2", col_lab = "QC2_cv")
