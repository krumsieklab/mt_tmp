# MetaboTools
#
# Correlate variable with dilution factors from quotient normalization.
# - for factors: will produce boxplot
# - for quantitative variable: will produce scatter plot
#
# last update: 2018-10-24
# authors: JK
#


mt_plots_qc_dilutionplot <- function(
  D,      # SummarizedExperiment input
  comp    # list of sample annotation column to compare with
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(length(comp)==1)
  
  # validate that there has been exactly one quotient normalization step
  q <- D%>% mti_res_get_path(c("pre","norm","quot"))
  if (length(q)>1) stop("There has been more than one quotient normalization call.")
  if (length(q)==0) stop("No quotient normalization performed.")
  # get dilution factors
  vd = q[[1]]$output$dilution
  
  # get variable to compare to
  if (!(comp %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", comp))
  vc = colData(D)[[comp]]
  
  # generate data frame for plotting
  dfplot <- data.frame(
    dilution.factor = vd,
    y = vc
  )
  colnames(dfplot)[2] <- comp
  
  # either produce boxplot or scatter plot
  if (is.character(vc) || is.factor(vc)) {
    # ensure factor
    dfplot[[comp]] = as.factor(dfplot[[comp]])
    # boxplot
    p <- dfplot %>% ggplot() +
      geom_boxplot(aes_string(x=comp,y="dilution.factor",color=comp)) + 
      ggtitle("quotient normalization dilution factors")
    
  } else {
    if (!is.numeric(vc)) stop(sprintf("'%s' has to be character, factor, or numeric.", comp))
    # scatter
    p <- dfplot %>% ggplot() +
      geom_point(aes_string(x=comp,y="dilution.factor")) + 
      ggtitle("quotient normalization dilution factors")
  }
  
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("dilution factor plot, '%s'",comp),
      output = list(p)
    )
  
  # return
  D
  
}

