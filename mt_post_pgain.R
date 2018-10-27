################################################################################
## MULTIPLE TESTING CORRECTION
################################################################################
#' mt_post_pgain
#'
#' add p-gain to metabolite-ratio tests
#'
#' @author Jonas Zierer
#' @import SummarizedExperiment
#' @importFrom dplyr %>% mutate gather left_join filter arrange
#' @param D SummarizedExperiment object
#' @param statname name of the statistical comparison
#' @param pcolumn name of p-value column
#' @return SummarizedExperiment modifies stats object
#' @export mt_post_pgain
mt_post_pgain <- function(D,
                             statname,
                             pcolumn = p.value,
                             ...){
    pcolumn <- enquo(pcolumn)

    ## are these ratio results?
    if (length(mti_res_get_path(D, c("modify","ratios"))) != 1)
        stop("must supply ratios to calculate p-gains")
    
    ## stat
    if(missing(statname))
        stop("statname must be given")

    ## find entry
    stat_id <- metadata(D)$results %>%
                         map_lgl(~"stats" %in% .x$fun && .x$output$name == statname) %>%
                         which()
    if(length(stat_id) == 0)
        stop("stat element with name ", statname, " does not exist")
    if(length(stat_id)  > 1)
        stop("there are multiple stat elements with name ", statname)

    ## DO CORRECTION
    rd <- rowData(D) %>%
        as.data.frame() %>%
        mutate(var = rownames(D)) %>%
        select(var, m1, m2)
    res <- metadata(D)$results[[ stat_id ]]$output$table %>%
                     left_join(rd, by = "var")
    res <- res %>%
        left_join(res %>% filter(is.na(m2)) %>% transmute(m1 = var, p_single_1 = p.value), by = "m1") %>%
        left_join(res %>% filter(is.na(m2)) %>% transmute(m2 = var, p_single_2 = p.value), by = "m2") %>%
        mutate(pgain = pmin(p_single_1, p_single_2) / p.value) %>%
        select(-m1, -m2, -p_single_1, -p_single_2)

    ## update results table
    metadata(D)$results[[ stat_id ]]$output$table <- res
    metadata(D)$results[[ stat_id ]]$output$logtxt %<>% str_c(., ", pgain")
        
    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Calculate p-gains for '%s'", statname),
                      output = NULL
                  )
    ## RETURN
    D
    
}
