#' Generate p-value histogram.
#'
#' Either for a given list of statistical results, or all.
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_names Name of the statistical results. If not given, will generate histogram for all.
#'
#' @return $result: plot, p-value histogram
#'
#' @examples
#' \dontrun{... %>% mt_plots_pvalhist() %>% ...                  # for all
#' ... %>% mt_plots_pvalhist(stat_names='comp') %>% ...  # for one
#' }
#'
#' @author JK
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_pvalhist <- function(
  D,
  stat_names=NULL
) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # if no stat_names given -> get all
  if (is.null(stat_names)) stat_names <- D %>% mti_res_get_stats_entries() %>% purrr::map("output") %>% purrr::map("name") %>% unlist()

  # loop over stat_names
  plots <- lapply(stat_names, function(statname){
    # breaks
    st <- D %>% mti_get_stat_by_name(statname)
    breaks <- pretty(range(st$p.value), n = grDevices::nclass.FD(st$p.value), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    # histogram
    p <- st %>%
        ggplot(aes(x=p.value)) +
        geom_histogram(binwidth=bwidth,fill="white",colour="black") +
        xlim(0,1) +
        ggtitle(glue::glue("'{statname}' p-values"))
    p
    })




  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("P-value histograms for {paste0(stat_names, collapse=', ')}"),
      output = plots
    )

  # return
  D

}


