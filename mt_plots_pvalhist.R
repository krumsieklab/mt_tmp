# MetaboTools
#
# P-value histogram(s)
# either for a given list of comparisons, or all
#
# last update: 2018-10-13
# authors: JK
#

require(glue)

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


