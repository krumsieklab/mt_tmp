### -- MT Reorder Samples Annotations for Plotting -- ###
# MetaboTools allows users to reorder sample annotations for plotting using the function mt_modify_sampanno_reorder()

library(MetaboTools)

D <-
  # load data
  mt_files_load_metabolon(file = codes.makepath("Mt/sampledata.xlsx"), sheet = "OrigScale") %>%
  # timing start
  mt_logging_tic() %>%

  ###
  # heading
  mt_reporting_heading("Sample Annotation Re-ordering") %>%
  mt_reporting_text("Users can reorder specific sample annotations using the function mt_modify_sampleanno_reorder. Annotations remain
                    in this new order for the remainder of the pipeline.") %>%
  # sample boxplot
  mt_plots_sampleboxplot(color=Group) %>%
  # filter metabolites with >20% missing values, then samples with >10% missing values
  mt_pre_filtermiss(met_max=0.2) %>%
  mt_pre_filtermiss(sample_max=0.1) %>%
  # batch correction by variable BATCH_MOCK
  mt_pre_batch_median(batches = "BATCH_MOCK") %>%
  # heading
  mt_reporting_heading("Part 2", lvl=2) %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  mt_plots_qc_dilutionplot(comp="Group") %>%
  # Reorder groups
  mt_modify_sampleanno_reorder(col_name = "Group", new_order = c("Vehicle", "treatment1", "treatment2")) %>%
  # dilution plot after reordering
  mt_plots_qc_dilutionplot(comp = "Group") %>%
  # final sample boxplot (with reodering)
  mt_plots_sampleboxplot(color=Group, plottitle = 'final')


D %>% mt_reporting_html(outfile="example_sample_reorder.html", output.calls = T)
