#' Compute p-gain from metabolite ratio test
#'
#' Add p-gain to result table
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name name of the statistical comparison
#' @param p_col name of p-value column to compute p-gain from
#'
#' @return $result: statistical object
#'
#' @examples
#' \dontrun{# add p-gains to the result table of the statistical comparison called "comparison1"
#' ... %>%
#'  mt_post_pgain(stat_name="comparison1") %>% ...
#'  }
#'
#' @author JZ
#'
#' @export
mt_post_pgain <- function(D,
                             stat_name,
                             p_col = p.value){
    p_col <- enquo(p_col)

    ## are these ratio results?
    if (length(mti_res_get_path(D, c("modify","ratios"))) != 1)
        stop("must supply ratios to calculate p-gains")

    ## stat
    if(missing(stat_name))
        stop("stat_name must be given")

    ## find entry
    stat_id <- metadata(D)$results %>%
                         map_lgl(~"stats" %in% .x$fun && .x$output$name == stat_name) %>%
                         which()
    if(length(stat_id) == 0)
        stop("stat element with name ", stat_name, " does not exist")
    if(length(stat_id)  > 1)
        stop("there are multiple stat elements with name ", stat_name)

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
                      logtxt = sprintf("Calculate p-gains for '%s'", stat_name),
                      output = NULL
                  )
    ## RETURN
    D

}
