#' mtw_analysis
#'
#' Creates an analysis result based on the type of analysis function chose by users, with plots for visulization
#' 
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name Name under which this comparison will be stored, must be unique to all other statistical results.
#' @param test_type Type of correlation to run, must be either "lm" or "cor". Default is "lm".
#' @param is_bin Logical. Whether the data is binary or not. Default is TRUE
#' @param lm_options A list of parameters for \code{mt_stats_univ_lm()}, expressions should be quoted by \code{quote()}
#' @param cor_options A list of parameters for \code{mt_stats_univ_cor()}, expressions should be quoted by \code{quote()}
#' @param multT_options A list of parameters for \code{mt_post_mulTest()}, expressions should be quoted by \code{quote()}
#' @param stat_metab_filter Filter formula for metabolites, used as "stat_filter" in \code{mt_logging_statsinfo()} and "metab_filter" in \code{mt_plots_volcano()}, should be quoted by \code{quote()}. Default is \code{quote(p.adj < 0.05)}
#' @param volcano_options A list of parameters for \code{mt_plots_volcano()}, expressions and variable names should be quoted by \code{quote()}
#' @param box_scatter_options A list of parameters for \code{mt_plots_boxplot_scatter()}, expressions and variable names should be quoted by \code{quote()}
#' 
#'
#' @return D with analysis results
#' @examples
#' \dontrun{ ... %>% mtw_analysis(
#' stat_name = 'Diag_mets',
#' test_type = 'cor',
#' is_bin = F,
#' cor_options = list(method = 'kendall', var = 'Diag_num', stat_name = 'Diag_mets'),
#' multT_options = list(method = 'BH'),
#' stat_metab_filter = quote(p.adj < 0.05)),
#' volcano_options = list(x = quote(statistic), colour = quote(p.adj < 0.05)),
#' box_scatter_options = list(plot_type = 'box',
#'                            x                  = quote(Diagnosis),
#'                            fill               = quote(Diagnosis),
#'                            metab_sort         = quote(p.value),
#'                            annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
#'                            cols               = 3)) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mtw_analysis <-
  function(D,
           stat_name,
           test_type = "lm",
           is_bin = T,
           lm_options = list(),
           cor_options = list(),
           multT_options = list(),
           stat_metab_filter = quote(p.adj < 0.05),
           volcano_options = list(),
           box_scatter_options = list()
           ) {
    stopifnot("SummarizedExperiment" %in% class(D))
    
    # help function to merge default and user input
    map_lists <- function(def, user) {
      # check if there are any illegal entries
      bad_args <- setdiff(names(user), names(def))
      if (length(bad_args) > 0) {
        stop(sprintf(
          "Invalid argument(s): %s",
          paste0(bad_args, collapse = ", ")
        ))
      }
      # now fill up entries in def
      def[names(user)] <- user
      # return
      def
    }
    
    lm_options_def = list(formula = '~ Age',stat_name = stat_name, sample_filter = NULL, n_cores = 1)
    
    cor_options_def = list(
      method = 'kendall',
      var = 'x',
      stat_name = stat_name,
      sample_filter = NULL,
      exact_flag = NULL
    )
    
    # Default parameters of mt_post_multTest
    # Removed "p_col = p.value"
    multT_options_def = list(stat_name = stat_name, method = "bonferroni")
    
    logging_stats_options = list(stat_name = stat_name, stat_filter = stat_metab_filter)
    
    volcano_options_def = list(
      stat_name = stat_name,
      x = quote(fc),
      metab_filter = stat_metab_filter,
      #xlabel = quote(gsub("~", "", as.character(x))),
      vline = NA,
      ggadd = NULL,
      colour = stat_metab_filter
    )
    
    box_scatter_options_def = list(
      plot_type = 'box',
      stat_name = stat_name,
      x = 'x',
      metab_filter = stat_metab_filter,
      metab_sort = quote(p.value),
      annotation = "{sprintf('P-value: %.1e', p.value)}",
      cols = 3,
      full_info = F,
      text_size = 3.88,
      jitter = "beeswarm",
      restrict_to_used_samples = T,
      manual_ylab = NULL,
      fitline = T,
      fitline_se = T,
      ggadd = NULL,
      fill = NULL
    )
    
    
    if (!(test_type == "lm")&&!(test_type == "cor")) {
      stop("Unknown test type")
    }
    
    # not able to make sample_filter as input param
    if (test_type == "lm") {
      #add statement, if sample filter is NULL or NA then delete it
      lm_options <- map_lists(lm_options_def, lm_options)
      if(is.null(lm_options$sample_filter)||is.na(lm_options$sample_filter)){lm_options$sample_filter <- NULL}
      lm_options$D <- D
      D <- do.call('mt_stats_univ_lm', lm_options)
    }
    if (test_type == "cor") {
      #add statement, if sample filter is NULL or NA then delete
      cor_options <- map_lists(cor_options_def, cor_options)
      if(is.null(cor_options$sample_filter)||is.na(cor_options$sample_filter)){cor_options$sample_filter <- NULL}
      cor_options$D <- D
      D <- do.call('mt_stats_univ_cor',cor_options)
    }
    
    
    # add fold changes to result tables
    if (is_bin == T) {
      D %<>% mt_post_addFC(stat_name = stat_name)
    }
    
    # add multiple testing correction
    multT_options <- map_lists(multT_options_def, multT_options)
    multT_options$D <- D
    D <- do.call('mt_post_multTest', multT_options)
        
    D %<>%
      # p-value histogram
      mt_plots_pvalhist(stat_names = stat_name) %>%
      mt_plots_pvalqq(stat_name = stat_name)
    
    
    # stat_name, stat_filter
    logging_stats_options$D <- D
    D <- do.call('mt_logging_statsinfo', logging_stats_options)
    
    
    #Volcano plot as overview of results
    volcano_options <- map_lists(volcano_options_def, volcano_options)
    volcano_options$D <- D
    D <- do.call('mt_plots_volcano', volcano_options)
    
    
    # (lm and binary) or cor -> boxplot; lm not binary -> scatter
    box_scatter_options <- map_lists(box_scatter_options_def, box_scatter_options)
    if(is.null(box_scatter_options$fill)||is.na(box_scatter_options$fill)){box_scatter_options$fill <- NULL}
    box_scatter_options$D <- D
    D <- do.call('mt_plots_boxplot_scatter', box_scatter_options)
    
  }