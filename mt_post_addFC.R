################################################################################
## ADD FOLDCHANGES TO RESULTS TABLE
################################################################################
#' mt_post_addFC
#'
#' add metabolite foldchanges to results table'
#' new column i called 'fc'
#'
#' @author Jonas Zierer
#' @import SummarizedExperiment
#' @param D SummarizedExperiment object
#' @param statname name of the statistical comparison
#' @param correct_confounder confounders to adjust for before calculating FC
#' @return SummarizedExperiment modifies stats object
#' @export mt_post_addFC
mt_post_addFC <- function(D,
                          statname,
                          correct_confounder,
                          combine = mean,
                          ...){

    ## FOLDCHANGE FUNCTION (CONSIDER PREVIOUS LOG)
    if (length(mti_res_get_path(D, c("pre","trans","log"))) != 1)
        stop("fold-changes can only be calculated for log-scale data")

    ## stat
    if(missing(statname))
        stop("statname must be given")
    ## FIND ENTRY
    stat_id <- metadata(D)$results %>%
                         map_lgl(~"stats" %in% .x$fun && .x$output$name == statname) %>%
                         which()
    if(length(stat_id) == 0)
        stop("stat element with name ", statname, " does not exist")
    if(length(stat_id)  > 1)
        stop("there are multiple stat elements with name ", statname)

    ## GET OUTCOME
    model <- metadata(D)$results[[ stat_id ]] %>%
                       .$output %>%
                       ## get first model in list (doesn't matter which)
                       .$lstobj %>%
                       .[[1]] %>%
                       ## get model matrix
                       .$model %>%
                       ## remove metabolites
                       select(-one_of(intersect(names(.), rownames(D)))) %>%
                       ## select first phenotype = outcome of interest
                       .[,1, drop = F]
    outcome    <- colnames(model)[[1]]

    ## CORRECT FOR CONFOUNDER
    formula    <- metadata(D)$results[[ stat_id ]]$output$formula
    confounder <- update.formula(formula, str_c(".~.-", outcome))
    confounder_terms <- attr(terms(confounder), "term.labels")
    if(length(confounder_terms) > 0){
        logmsg(glue::glue("correcting for {confounder}"))
        D <- mti_correctConfounder(D, confounder)
    }

    ## EXTRACT DATA
    d_fc <- mti_format_se_samplewise(D) %>%
        ## use only levels that were used in stat
        inner_join(model, by = outcome)

    ## CHECK TYPE OF OUTCOME
    if(!(class(d_fc[[ outcome ]]) %in% c("factor", "character")))
        stop(glue::glue("Fold-changes are only meaningful for factor/character, but {outcome} is a {class(d[[outcome]])}"))
    if(length(unique(d_fc[[ outcome ]])) != 2)
        stop(glue::glue("Fold-changes are only meaningful 2 groups, but {outcome} has {length(unique(d[[outcome]]))}"))

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
        gather(var, value, one_of(rownames(D))) %>%
        group_by(var,!!outcome_q) %>%
        summarise(value = combine(value)) %>%
        ungroup() %>%
        spread(!!outcome_q, value) %>%
        mutate(fc = !!levels_1 - !!levels_2) %>%
        select(var, fc)

    ## ADD TO RESULTS
    metadata(D)$results[[stat_id]]$output$table %<>%
                  left_join(d_fc, by = "var")

    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Calculate foldchanges for %s", statname),
                      output = NULL
                  )
    ## RETURN
    D
}

