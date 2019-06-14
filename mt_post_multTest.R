
#' Multiple testing correction
#' 
#' Adjust output of statistical tests for multiple testing
#'
#' @param D \code{SummarizedExperiment} input
#' @param statname name of the statistical comparison to adjust
#' @param pcolumn name of p-value column to adjust
#' @param method which method to use for multiple testing (see p.adjust)
#' 
#' @return $result: statistical object
#' 
#' @examples
#' # correct the statistical comparison called "Li's" using Benjamini-Hochberg
#' ... %>%
#'  mt_post_multTest(statname="Li's", method="BH") %>% ...
#' 
#' @author JK, JZ
#' 
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
                         map_lgl(~"stats" %in% .x$fun && .x$output$name == statname) %>%
                         which()
    if(length(stat_id) == 0)
        stop("stat element with name ", statname, " does not exist")
    if(length(stat_id)  > 1)
        stop("there are multiple stat elements with name ", statname)

    ## DO CORRECTION
    metadata(D)$results[[stat_id]]$output$table %<>%
                  mutate(p.adj = p.adjust(!!pcolumn, method = method))
    ## UPDATE LOG ENTRIES
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
