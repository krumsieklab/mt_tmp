
#' mt_plot_volcano
#'
#' create volcano plot
#'
#' @author Jonas Zierer, Jan Krumsiek
#' @import SummrizedExperiment
#' @import ggplot2
#' @importFrom dplyr %>% mutate gather left_join filter arrange
#' @importFrom rlang enquo !!
#' @param D SummarizedExperiment object
#' @param x what value shall be plotted on x (default: fc)
#' @param statname name of the statistics obkect to plot
#' @param metab_filter if given, filter will be applied to data and remaining varaibles will be labelled in plot
#' @param ... further parameters forwarded to ggplot::aes
#' @return SummarizedExperiment with volcano plot in metadata(D)$results
#' @export mt_plot_volcano
mt_plots_volcano <- function(D,
                             x = fc,
                             statname,
                             metab_filter = p.value < 0.05,
                             xlabel=gsub("~","",as.character(x)),
                             ...){
    x <- enquo(x)
    
    ## check input
    stopifnot("SummarizedExperiment" %in% class(D))
    if(missing(statname))
        stop("statname must be given for volcanoplot")

    ## rowData
    rd <- rowData(D) %>%
        as.data.frame() %>%
        mutate(var = rownames(D))
    
    ## stat
    data_plot <- mti_get_stat_by_name(D, statname) %>%
        inner_join(rd, by = "var") %>%
        mutate(xxx = !!x)

    ## SCALE -log10
    reverselog_trans <- function (base = exp(1)){
        trans <- function(x) -log(x, base)
        inv <- function(x) base^(-x)
        scales::trans_new(paste0("reverselog-", format(base)), trans, inv, 
                          scales::log_breaks(base = base),
                          domain = c(1e-100, Inf))
    }
    
    
    ## CREATE PLOT
    p <- data_plot %>%
        ## do plot
        ggplot(aes(x = xxx, y = p.value)) +
        geom_point(aes(...)) +
        scale_y_continuous(trans = reverselog_trans(10),
                           breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        labs(x = xlabel, y = "p-value") +
        ggtitle(statname)

    ## ADD METABOLITE LABELS
    if(!missing(metab_filter)){
        mti_logstatus("add label")
        metab_filter_q <- enquo(metab_filter)
        data_annotate <- data_plot %>%
            filter(!!metab_filter_q)
        p <- p + ggrepel::geom_text_repel(data = data_annotate,
                           aes(label = name))
    }
    
    ## ADD AXIS GROUPS
    d <- mti_get_stat_by_name(D, statname, fullstruct=T)
    if ("groups" %in% names(d) && length(d$groups)==2) {
      p <- mti_add_leftright_gg(p, paste0(d$groups[1],' high'), paste0(d$groups[2],' high'))
    }

    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("volcano plot, aes: %s", mti_dots_to_str(...)),
                      output = list(p)
                  )
    ## return
    D
}
