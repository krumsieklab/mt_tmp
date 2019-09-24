require(formula.tools)

#' Univariate GLMs.
#' 
#' Computes univariate GLM for each metabolite.
#' 
#' \enumerate{
#'   \item Will treat the first term of the formula as outcome.
#'   \item If outcome has >2 factors, will perform ANOVA.
#'   \item If random effect term is contained, will used lmer function. (can also be used for paired analysis, e.g. before/after)
#' }
#' 
#' @param D \code{SummarizedExperiment} input
#' @param formula left-hand side of formula to be put into glm function.
#' @param name name under which this comparison will be stored, must be unique to all other statistical results
#' @param samplefilter term which samples to filter to first... e.g. used if the data contains >2 groups but the user wants to run a two-group comparison
#' @param mc.cores number of cores to use for mclapply... default: 1. More than one core will not work on Windows platforms.
#'
#' @return $result: statistics object
#'
#' @examples
#' # run lm with no confounders, "Group" as outcome
#' # filter to groups "Li_2" and "Li_5"
#' # name the comparison "Li's"
#' ... %>%
#'  mt_stats_univ_lm(
#'    formula      = ~ Group, 
#'    samplefilter = (Group %in% c("Li_2","Li_5")),
#'    name         = "Li's"
#'  ) %>% ...
#' 
#' @author JK, JZ
#' 
mt_stats_univ_lm <- function(
  D,              # SummarizedExperiment input
  formula,        # formula defining statistical model, see above.
  name,           # name of comparison,
  samplefilter,
  mc.cores = 1
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # make sure name does not exist yet
  if (name %in% unlist(mti_res_get_stats_entries(D) %>% map("output") %>% map("name"))) stop(sprintf("stat element with name '%s' already exists",name))
  
  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise() 
  
  ## FILTER SAMPLES
  if(!missing(samplefilter)) {
    
    filter_q <- enquo(samplefilter)
    Ds <- Ds %>%
      filter(!!filter_q) %>%
      droplevels()
    # message("filter metabolites: ", metab_filter_q, " [", nrow(stat), " remaining]")
    # did we leave 0 rows?
    if (nrow(Ds)==0) stop("Filtering left 0 rows")
    if (nrow(Ds)==ncol(D)) mti_logwarning('filtering did not filter out any samples')
    
  }
  
  ## save outcome variable
  outvar       <- attr(terms(formula, keep.order = T),"term.labels")[1]
  outvar_label <- outvar
  is_interaction <- FALSE
  if(str_detect(outvar, ":")){
    outvar <- strsplit(outvar, ":") %>% unlist()
    is_interaction <- TRUE
    mti_logmsg("calculating interaction term")
  }
  if(any(!(outvar %in% colnames(Ds))))
    stop(sprintf("column %s do not exist in data", str_c(outvar[ !(outvar %in% colnames(Ds))], collapse = ", ")))
  
  # handle factor outcomes
  do_anova    <- FALSE
  outvar_term <- outvar
  for(o in seq_along(outvar)){  
    v = outvar[o]
    if (is.character(Ds[[v]]))
      Ds[[ v ]] = as.factor(Ds[[ v ]])
    if (is.factor(Ds[[v]])) {
      # check that there are exactly two levels
      if (length(levels( Ds[[ v]] ))!=2){
        stop(sprintf("factor outcomes must have exactly two levels, '%s' has %d", o, length(levels(v))))
        do_anova <- TRUE
      }
      # remember the level name of the second (will be deleted later on)
      outvar_term[o] <- str_c(outvar[o], levels(Ds[[ v ]])[2])
    }
  }
  if(is_interaction)
    outvar_term <- str_c(outvar_term, collapse = ":")
  
  
  ## choose lm functions
  has_random_eff <- FALSE
  if(str_detect(as.character(formula), "\\|")){
    mti_logmsg("random effect detected; using lmer")
    f_lm     <- lmerTest::lmer
    f_tidy   <- broom.mixed::tidy
    has_random_eff <- TRUE
  }else{
    f_lm   <- lm
    f_tidy <- broom::tidy
  }
  f_tidy_tidy <- function(m, ...){
    if(is.null(m))
      return(tibble(term = outvar_term))
    f_tidy(m) %>%
      mutate(formula = as.character(m$terms)) %>%
      dplyr::select(term, formula, everything())
  }
  
  
  ## validate formula
  if (!is.null(lhs(formula))) stop("Left-hand side of formula must be empty")
  ## will crash sort of meaningfully if variables don't exist
  if(has_random_eff){
    mm <- lFormula(update.formula(formula, str_c(rownames(D)[[1]], "~.")), Ds)
  }else{
    mm <- model.matrix(formula,Ds)
  }
  
  ## check if anova is necessary for interaction effects
  ## (happens if 2 factors with 2 levels each where individual variables are not included)
  if(is_interaction){
    interactions <- map(outvar, function(.x){
      if(is.factor(Ds[[.x]]))
        str_c(.x, unique(Ds[[.x]]))
      else
        .x
    })%>%
      expand.grid() %>%
      apply(MARGIN = 1, str_c, collapse = ":")
    if(length(intersect(colnames(mm), interactions)) > 1){
      do_anova <- TRUE
      stop("for interactions effects between factors, individual terms should be included")
    }
  }
  
  
  ## function to create models
  ## does checks to minimise errors  
  do_lm <- function(m){
    ## run glm with updated formula
    form <- update(formula, sprintf("%s~.",m))
    ## check for constant confounders
    trms <- attr(terms(formula), "term.labels") %>%
      purrr::discard(~str_detect(.x, ":")) %>%
      c(m)
    clss <- map_chr(trms, ~class(Ds[[.x]])) %>%
      setNames(trms)
    ## subset to complete data
    d <- Ds %>%
      dplyr::select(one_of(trms), !!rlang::sym(m)) %>%
      filter(complete.cases(.))
    ## check for invariant vonfounders
    conf_invar_num <- clss %>%
      purrr::keep(~.x %in% c("integer", "numeric")) %>%
      imap(~var(d[[.y]])) %>%
      purrr::keep(~.x == 0)
    conf_invar_fct <- clss %>%
      purrr::discard(~.x %in% c("integer", "numeric")) %>%
      imap(~length(unique(d[[.y]]))) %>%
      purrr::keep(~.x == 1)
    conf_invar <- c(conf_invar_num, conf_invar_fct)
    if(length(conf_invar) > 0){
      ## terminate if metabolite or outcome are invariant
      if(m %in% names(conf_invar)){ 
        mti_logwarning(glue::glue("metabolite {m} invariant "))
        return(NULL)
      }
      if(any(outvar %in% names(conf_invar))){ 
        mti_logwarning(glue::glue("outcome {outvar} invariant for metabolite {m}"))
        return(NULL)
      }
      for(c in names(conf_invar)){
        mti_logwarning(glue::glue("confounder {c} invariant for metabolite {m}, removing from formula"))
        form <- update.formula(form, str_c(". ~ . -", c))
      }
    }
   
    ## CALCLUATE ACTUAL MODEL
    mod <- f_lm(
      data    = Ds,
      formula = form
    )
    ## DO ANOVA IF MULTIPLE FACTOR LEVELS
    if(do_anova)
      mod <- anova(mod)
    ## RETURN
    mod
  }
  
  ## run tests for all metabolites
  models <- parallel::mclapply(rownames(D), do_lm, mc.cores = mc.cores) %>%
    setNames(rownames(D))
  
  # broom it up, subselect to term, rename term
  tab <- map_dfr(models, f_tidy_tidy, conf.int = T, .id = "var") %>%
    dplyr::filter(term == outvar_term) %>%
    mutate(term =  outvar_label)
  
  ## tidy up a bit more
  tab <- tab %>%
    dplyr::select(-matches("^(effect|group)$"))
  
  ## construct output groups variable
  if (is.factor(Ds[[outvar]])) {
    outgroups <- levels(Ds[[outvar]])
  } else {
    outgroups <- NULL
  }
  
  ## add status information & results
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("univariate lm, %s", as.character(formula)),
      output = list(
        table   = tab,
        formula = formula,
        name    = name,
        lstobj  = models,
        groups = outgroups
      )
    )
  
  ## return
  D
}
