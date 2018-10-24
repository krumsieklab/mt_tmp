# MetaboTools
#
# Univariate GLMs, one per metabolite.
#
# model: met ~ outcome + conf1 + conf2
# input formula as this: ~ outcome + conf1 + conf2
# first variable on right-hand side of formula is always considered as the outcome.
#
# Note on outcome variable:
# - cannot handle interaction term outcomes right now
# - factor outcomes can only have exactly 2 levels
#
# last update: 2018-10-23
# authors: JK
#

require(formula.tools)

mt_stats_univ_lm <- function(
  D,              # SummarizedExperiment input
  formula,        # formula defining statistical model, see above.
  name = '',  # [optional] name of comparison,
  samplefilter
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise() 
  
  ## FILTER METABOLITES
  if(!missing(samplefilter)) {
    filter_q <- enquo(samplefilter)
    Ds <- Ds %>%
      filter(!!filter_q) %>% droplevels()
    # message("filter metabolites: ", metab_filter_q, " [", nrow(stat), " remaining]")
    # did we leave 0 rows?
    if (nrow(Ds)==0) stop("Filtering left 0 rows")
    if (nrow(Ds)==ncol(D)) warning('filtering did not filter out any samples')
  }

  # save outcome variable
  outvar <- attr(terms(formula),"term.labels")[1]
  if(!(outvar %in% colnames(Ds))) stop(sprintf("column %s does not exist in data",outvar))
  
  # handle factor outcomes
  v = Ds[[outvar]]
  termsuff = ''
  if (is.character(v) || is.factor(v)) {
    # convert from string
    if (is.character(v)) v=as.factor(v)
    # check that there are exactly two levels
    if (length(levels(v))!=2) stop(sprintf("factor outcomes must have exactly two levels, '%s' has %d", outvar, length(levels(v))))
    # remember the level name of the second (will be deleted later on)
    termsuff = levels(v)[2]
  }

  # validate formula
  if (!is.null(lhs(formula))) stop("Left-hand side of formula must be empty")
  model.matrix(formula,Ds) # will crash sort of meaningfully if variables don't exist
  
  
  # iterate over data matrix, run tests
  X = t(assay(D))
  models <- lapply (rownames(D), function(f) {
     # run glm with updated formula
    lm(
      data = Ds,
      formula = update(formula, sprintf("%s~.",f))
      )
  })
  
  # broom it up, subselect to term, rearrange, possibly rename factor field
  tab <- map_dfr(models, broom::tidy) %>% filter(term==paste0(outvar,termsuff)) %>% 
     mutate(var=rownames(D)) %>% select(var,everything()) %>% select(-term,term)
  tab$term = outvar

  # set attributes
  attr(tab, 'formula') <- formula
  attr(tab, 'name')    <- name
  attr(tab, "lstobj")  <- models
  
  
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


