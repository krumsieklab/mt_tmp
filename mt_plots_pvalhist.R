require(glue)
#' Generate p-value histogram.
#' 
#' Either for a given list of statistical results, or all.
#'
#' @param D \code{SummarizedExperiment} input
#' @param statnames Name of the statistical results. If not given, will generate histogram for all.
#'
#' @return $result: plot, p-value histogram
#'
#' @examples
#' ... %>% mt_plots_pvalhist() %>% ...                  # for all
#' ... %>% mt_plots_pvalhist(statnames='comp') %>% ...  # for one
#' 
#' @author JK
mt_plots_pvalhist <- function(
  D,
  statnames=NULL
) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # if no statnames given -> get all
  if (is.null(statnames)) statnames <- D %>% mti_res_get_stats_entries() %>% map("output") %>% map("name") %>% unlist()
  
  # loop over statnames
  plots <- lapply(statnames, function(statname){
    # breaks
    st <- D %>% mti_get_stat_by_name(statname) 
    breaks <- pretty(range(st$p.value), n = nclass.FD(st$p.value), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    # histogram
    st %>% 
        ggplot(aes(x=p.value)) +
        geom_histogram(binwidth=bwidth,fill="white",colour="black") +
        xlim(0,1) +
        ggtitle(glue("'{statname}' p-values"))
    })
  

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue("P-value histograms for {paste0(statnames, collapse=', ')}"),
      output = plots
    )
  
  # return
  D
  
}


