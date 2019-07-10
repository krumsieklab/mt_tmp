#### In here we collect examples how to access the internal structures of MT, primarily the metadata(D)$output structure.



#### load MT ----
mt.quickload()

#### part 1, run some standard pipeline, no questions asked ----

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
  mt_pre_impute_knn() %>%
  # scale %>% 
  mt_pre_trans_scale() %>% 
  # PCA
  mt_plots_PCA(labelby = 'CLIENT_IDENTIFIER', color='GROUP_DESC') %>%
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group, 
    samplefilter = (Group %in% c("treatment1","treatment2")),
    name         = "comp",
    mc.cores     = 1
  ) %>%
  # linear model, differential test on Group
  mt_stats_univ_lm(
    formula      = ~ Group, 
    samplefilter = (Group %in% c("Vehicle","treatment2")),
    name         = "comp2",
    mc.cores     = 1
  ) 


#### part 2, examples for internal access ----

## example: access dilution factor vector

# 'non-elegant' way, access via index
metadata(D)$results[[11]]
metadata(D)$results[[11]]$output

# more elegant way, access via name
# careful: could return multiple results
L <- D %>% mti_res_get_path(c('pre','norm','quot'))
# only access the first (in this case, there only is one)
L[[1]]


## example: get all dilution plots
L <- D %>% mti_res_get_path(c('plots','qc','dilutionplot'))
# access one, then the other, change to BW background
L[[1]]$output[[1]] + theme_bw()
L[[2]]$output[[1]] + theme_bw()


## example: access data frame with statistical results
res <- D %>% mti_get_stat_by_name("comp2")

# show
res
# add nice names from rowData
rd <- rowData(D) %>%
  as.data.frame() %>%
  mutate(var = rownames(D))
res %>% inner_join(rd, by = "var") # joint data frame of results and rowData


