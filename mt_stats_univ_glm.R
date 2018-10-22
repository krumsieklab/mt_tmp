# MetaboTools
#
# Univariate GLMs, one per metabolite.
#
#
# last update: 2018-10-21
# authors: JK
#
# WARNING
# - currently only does lm
# - still not fully safe implementation (handling of formula, subsetting to 2nd result entry after brooming)

mt_stats_univ_glm <- function(
  D,              # SummarizedExperiment input
  formula,        # formula defining statistical model
  compname = ''   # [optional] name of comparison
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # replace first argument in formula with %s placeholder
  dp <- deparse(formula)
  if (length(grep("\\+", dp))>0) {
    sform <- sub("~.*? \\+", "~ %s +", dp)
  } else {
    sform <- sub("~.*", "~ %s", dp)
  }
  
  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise() 
  
  # iterate over data matrix
  X = t(assay(D))
  models <- lapply (1:ncol(X), function(i) {
    # run glm
    glm(
      data = Ds,
      formula = sprintf(sform, colnames(X)[i])
      )
  })
  # get all second lines (stats of our metabolite)
  tab <- map_dfr(models, function(x){broom::tidy(x)[2,]})
  colnames(tab)[colnames(tab)=="term"] <- "var"

  # set attributes
  attr(tab, 'formula') <- formula
  attr(tab, 'name') <- compname
  attr(tab, "lstobj") <- models
  
  
  # add status information & results
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("univariate lm, %s",deparse(formula)),
      output = tab
    )
  
  # return
  D
  
  
  
}




