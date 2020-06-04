
#### load MT ----
library(MetaboTools)

#### part 1, needed preprocessing steps ----

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
  mt_pre_filtermiss(met_max=0.2) %>%
  mt_pre_filtermiss(sample_max=0.1) %>%
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
  mt_pre_trans_scale()


#### step 2, new functions ----

D1 <- D %>%
  mt_files_write_xls(file = "out.xlsx") %>%
  mt_pre_outlier(method="univariate", thresh=1, perc=0.12) %>%
  # mt_plots_PCA() %>%
  mt_plots_PCA(labelby = 'CLIENT_IDENTIFIER', color='GROUP_DESC', expvarplot = T) %>%
  mt_plots_PCA(labelby = 'CLIENT_IDENTIFIER', color='outlier_univariate') %>%
  mt_plots_PCA(labelby = 'CLIENT_IDENTIFIER', color='GROUP_DESC', expvarplot = T) %>%
  mt_plots_PCA(labelby = 'CLIENT_IDENTIFIER', color='GROUP_DESC', ellipse = 0.95) %>%
  mt_plots_PCA(PCa = 1, PCb = 3) %>%
  mt_plots_PCA(show = 'loadings') %>%
  mt_plots_PCA(show = 'loadings', PCa = 1, PCb = 3) %>%
  mt_plots_PCA(show = 'loadings', labelby = 'name', textrepel = F) %>%
  mt_plots_PCA(show = 'loadings', labelby = 'name', textrepel = F, color='SUPER_PATHWAY') %>%
  # filter out outliers, and redo PCA
  mt_modify_filter_samples(sample_filter = outlier_univariate == FALSE) %>%
  mt_plots_PCA()


#### step 3, create output ----
D1 %>% mt_reporting_html("example3.html")
# D1 %>% mt_reporting_generateMD(outfile = "example3.RMD")

# # mutate examples, stored
# D %>% mt_modify_mutate(dir='samples', varname='new', term=num2^2) %>% colData()
D %>% mt_modify_mutate(anno_type='samples', col_name='new', term= (grepl("Vehicle", Group))) %>% colData()

# D %>% mt_modify_mutate(dir='metabolites', varname='new2', term=MASS-100) %>% rowData()

... %>% mt_modify_mutate(dir='samples', varname='is_studypool', term=(grepl("study", ID))) %>% ...
