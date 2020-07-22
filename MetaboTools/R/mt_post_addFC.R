#' Compute fold-change
#'
#' Add metabolite fold-changes to results table
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name name of the statistical comparison
#' @param correct_confounder confounders to adjust for before calculating fold-changes
#' @param combine function with respect to which fold changes are computed (e.g. difference of the log means or median)
#'
#' @return $result: statistical object
#'
#' @examples
#' \dontrun{# add fold-changes to the result table of the statistical comparison called "comparison1", after correcting for variable "age"
#' ... %>%
#'  mt_post_addFC(stat_name="comparison1", correct_confounder=age) %>% ...
#'  }
#'
#'
#' @author JZ
#'
#' @export
mt_post_addFC <- function(D,
                          stat_name,
                          correct_confounder,
                          combine = function(x){mean(x,na.rm=T)}){

  ## FOLDCHANGE FUNCTION (CONSIDER PREVIOUS LOG)
  if ((length(MetaboTools:::mti_res_get_path(D, c("pre","trans","log"))) != 1) &&
      (length(MetaboTools:::mti_res_get_path(D, c("flag","logged"))) != 1))
    stop("fold-changes can only be calculated for log-scale data")

  ## stat
  if(missing(stat_name))
    stop("stat_name must be given")
  ## FIND ENTRY
  stat_id <- metadata(D)$results %>%
    purrr::map_lgl(~"stats" %in% .x$fun && .x$output$name == stat_name) %>%
    which()
  if(length(stat_id) == 0)
    stop("stat element with name ", stat_name, " does not exist")
  if(length(stat_id)  > 1)
    stop("there are multiple stat elements with name ", stat_name)

  ## stop if the results are coming from an ANOVA, determine from first model
  if (class( metadata(D)$results[[ stat_id ]] %>%
             .$output %>%
             ## get first model in list (doesn't matter which)
             .$lstobj %>%
             .[[1]])[[1]] == "anova") {
    stop("Cannot add fold change for ANOVA analysis.")
  }

  ## GET OUTCOME (works for both lm and lmer)
  formula <- metadata(D)$results[[ stat_id ]]$output$formula
  terms <- all.vars(as.formula(formula))
  outcome <- terms[1]

  # find any random effects... they cannot be corrected for here
  # parse out by finding any variable inside (x | y)
  re_terms <- stringr::str_match(formula, "\\((.*?)\\|(.*?)\\)")[-1] %>% trimws()
  terms <- setdiff(terms, re_terms)


  ## CORRECT FOR CONFOUNDER
  if (length(terms)>1) {
    # construct new formula
    conf_form <- as.formula(sprintf("~%s", paste0(terms[2:length(terms)], collapse = "+")))
    MetaboTools:::mti_logstatus(glue::glue("correcting for {as.character(conf_form)}"))
    D <- MetaboTools:::mti_correctConfounder(D, conf_form)
  }

  ## EXTRACT DATA
  # if(is.factor(colData(D)[[outcome]]))model[[outcome]] <- as.factor(model[[outcome]])
  d_fc <- cbind(colData(D)[,outcome,drop=F],t(assay(D))) %>% as.data.frame()
  keep <- metadata(D)$results[[ stat_id ]]$output$samples.used
  d_fc <- d_fc[keep,,drop=F]
  d_fc[[outcome]] <- as.factor(d_fc[[outcome]])
  d_fc[[outcome]] <- droplevels(d_fc[[outcome]])

  ## CHECK TYPE OF OUTCOME
  if(!(class(d_fc[[ outcome ]]) %in% c("factor", "character")))
    stop(glue::glue("Fold-changes are only meaningful for factor/character, but {outcome} is a {class(d_fc[[outcome]])}"))
  if(length(unique(d_fc[[ outcome ]])) != 2)
    stop(glue::glue("Fold-changes are only meaningful for 2 groups, but {outcome} has {length(unique(d_fc[[outcome]]))}"))

  ## GET LEVELS IN ORDER OF FACTOR LEVELS
  outcome_q <- rlang::sym(outcome)
  if( "factor" %in% class(d_fc[[outcome]])){
    levels <- levels(d_fc[[outcome]])
  }else{
    levels <- unique(d_fc[[outcome]])
  }
  levels_1  <- rlang::sym(levels[1])
  levels_2  <- rlang::sym(levels[2])

  ## FOLDCHANGE FUNCTION (CONSIDER PREVIOUS LOG)
  d_fc <- d_fc %>%
    tidyr::gather(var, value, one_of(rownames(D))) %>%
    dplyr::group_by(var,!!outcome_q) %>%
    dplyr::summarise(value = combine(value)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(!!outcome_q, value) %>%
    dplyr::mutate(fc = !!levels_1 - !!levels_2) %>%
    dplyr::select(var, fc)

  ## ADD TO RESULTS
  metadata(D)$results[[stat_id]]$output$table %<>%
    dplyr::left_join(d_fc, by = "var")

  ## make sure fold change has the same sign as statistic
  ## this is a debug solution, it should come out properly from the code above... but this fixes the bug for now
  metadata(D)$results[[stat_id]]$output$table$fc <-
    metadata(D)$results[[stat_id]]$output$table$fc %>% abs() * sign(metadata(D)$results[[stat_id]]$output$table$statistic)


  ## add status information & plot
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Calculated foldchanges for %s", stat_name),
      output = NULL
    )
  ## RETURN
  D
}

