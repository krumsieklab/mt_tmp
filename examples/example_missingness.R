#####
##### MetaboTools pipeline
#####

#### Missingness analysis in two groups

#### run pipeline ----

# load MT
library(MetaboTools)

D <-
  # load data
  mt_files_load_metabolon(codes.makepath("Mt/sampledata.xlsx"), "OrigScale") %>%
  # timing start
  mt_logging_tic() %>%

  # missingness analysis between treatments
  mt_stats_univ_missingness(
    comp_col = 'Group',
    stat_name = 'miss',
    sample_filter = (Group %in% c("treatment1","treatment2"))) %>%

  # add multiple testing correction
  mt_post_multTest(stat_name = "miss", method = "BH") %>%

  # Volcano plot (with statistic)
  mt_plots_volcano(stat_name     = "miss",
                   x = statistic,
                   metab_filter = p.value < 0.05,
                   colour       = p.value < 0.05) %>%

  # final timing
  mt_logging_toc() %>%

  # END OF PIPELINE
  {.}




#### generate and knit markdown ----

D %>% mt_reporting_html(outfile="example_missingness.html")
