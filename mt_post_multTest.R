################################################################################
## MULTIPLE TESTING CORRECTION
################################################################################
#' mt_post_multTest 
#'
#' adjust output of stats for multiple testing
#'
#' @author Jonas Zierer
#' @import SummarizedExperiment
#' @importFrom dplyr %>% mutate gather left_join filter arrange
#' @param D SummarizedExperiment object
#' @param statname name of the statistical comparison
#' @param method which method to use for multiple testing (see p.adjust)
#' @return SummarizedExperiment modifies stats object
#' @export mt_post_multTest
mt_post_multTest <- function(D,
                             statname,
                             pcolumn = p.value,
                             method = "bonferroni",
                             ...){
    pcolumn <- enquo(pcolumn)
    
    ## stat
    if(missing(statname))
        stop("statname must be given")

    ## FIND ENTRY
    stat_id <- metadata(D)$results %>%
                         map_lgl(~"stats" %in% .x$fun && attr(.x$output, "name") == statname) %>%
                         which()
    if(length(stat_id) == 0)
        stop("stat element with name ", statname, " does not exist")
    if(length(stat_id)  > 1)
        stop("there are multiple stat elements with name ", statname)

    ## DO CORRECTION
    metadata(D)$results[[stat_id]]$output %<>%
                  mutate(p.adj = p.adjust(!!pcolumn, method = method))
    metadata(D)$results[[stat_id]]$logtxt %<>%
                  str_c(., ", multiple testing adjusted (", method, ")")
    metadata(D)$results[[stat_id]]$logshort %<>%
                  str_c(., " [", method, "]")

    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Multiple testing correction of '%s' using '%s'", statname, method),
                      output = NULL
                  )
    ## RETURN
    D
    
}
