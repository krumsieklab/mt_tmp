#' mtw_apply_analysis
#'
#' Running standard analyses for multiple phenotypes, analysis steps including lm/cor analysis, 
#' multi-test with histgram, QQ-plot and volcano plot showing p-values, and box/scatter plot for the 
#' visualization of significant metabolites.
#' 
#'
#' @param D \code{SummarizedExperiment} input
#' @param outcomes A data frame, should contain two columns: outcomes$outcome: name of phenotype columns in colData; outcomes$type: data type of pheotype columns, must be one of numeric / binary / ordinal
#' 
#'
#' @return D with multiple analysis results and plots
#' @examples
#' \dontrun{... %>% mtw_apply_analysis(outcomes) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mtw_apply_analysis <- function(D, outcomes){
  
  #nrow(outcomes)
  for (i in 1:nrow(outcomes)) {
    outcome <- outcomes$outcome[i]
    outcome_type <- outcomes$type[i]
    outcome_name <- paste0(outcome, "_", outcome_type)
    
    input_list <- list(
      stat_name = quote(outcome_name),
      is_bin = F,
      multT_options = list(method = 'BH'),
      logging_stats_options = list(stat_filter = quote(p.adj < 0.05)),
      volcano_options = list(x = quote(statistic), metab_filter = quote(p.adj < 0.05), 
                             colour = quote(p.adj < 0.05))
    )
    
    # ordinal variables are actually modeled by their numeric counterparts and Kendall's tau
    if (outcome_type == 'ordinal') {
      D %<>% mt_reporting_heading(strtitle = paste0(outcome,'\n'), lvl = 2) %>% 
        mt_modify_mutate(anno_type = "samples", col_name = outcome, 
                         term = as.numeric(as.matrix(!!sym(outcome))))
      
        input_list_add <- list(
          D = D,
          test_type = 'cor',
          cor_options = list(method = 'kendall', var = outcome, stat_name = outcome_name),
          box_scatter_options = list(plot_type = 'box',
                                     x                  = sym(outcome),
                                     fill               = sym(outcome),
                                     metab_filter       = quote(p.adj < 0.05),
                                     metab_sort         = quote(p.value),
                                     annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                                     cols               = 3)
        )
      # analysis wrapper
      input_list[names(input_list_add)] <- input_list_add
      D <- do.call('mtw_analysis', input_list) # for numeric, twofactor, or multifactor: use lm
      
    } else if (outcome_type == 'numeric') {
      
      D %<>% mt_reporting_heading(strtitle = paste0(outcome,'\n'), lvl = 2) %>% 
        mt_modify_mutate(anno_type = "samples", col_name = outcome, 
                         term = as.numeric(as.matrix(!!sym(outcome))))
      
      input_list_add <- list(
        D = D,
        test_type = 'lm',
        lm_options = list(formula = formula(sprintf("~ %s",outcome))),
        box_scatter_options = list(plot_type = 'scatter',
                                   x                  = sym(outcome),
                                   fill               = NULL,
                                   metab_filter       = quote(p.adj < 0.05),
                                   metab_sort         = quote(p.value),
                                   annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                                   cols               = 3)
      )
      
      input_list[names(input_list_add)] <- input_list_add
      # analysis wrapper
      D <- do.call('mtw_analysis', input_list)
      
    } else if (outcome_type == 'binary'){
      
      D %<>% mt_reporting_heading(strtitle = paste0(outcome,'\n'), lvl = 2) %>% 
        mt_modify_mutate(anno_type = "samples", col_name = outcome, 
                         term = as.factor(!!sym(outcome)))
      
      input_list_add <- list(
        D = D,
        test_type = 'lm',
        is_bin = T,
        lm_options = list(formula = formula(sprintf("~ %s",outcome))),
        box_scatter_options = list(plot_type = 'box',
                                   x                  = sym(outcome),
                                   fill               = sym(outcome),
                                   metab_filter       = quote(p.adj < 0.05),
                                   metab_sort         = quote(p.value),
                                   annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}",
                                   cols               = 3)
      )
      
      # now fill up entries in def
      input_list[names(input_list_add)] <- input_list_add
      # analysis wrapper
      D <- do.call('mtw_analysis', input_list)
    }
  }
  return(D)
}