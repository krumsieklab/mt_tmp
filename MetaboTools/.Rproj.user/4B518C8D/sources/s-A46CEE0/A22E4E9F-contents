#' mt_plots_volcano
#'
#' Creates a volcano plot
#'
#' @param D \code{SummarizedExperiment} input
#' @param x what value shall be plotted on x (default: fc)
#' @param statname name of the statistics object to plot
#' @param metab_filter if given, filter will be applied to data and remaining varaibles will be labelled in plot
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param vline where to draw vertical line (for fold-change), has to be single value. default: none
#' @param hline where to draw horizontal line (for p-values), has to be an expression such as 'p.adj < 0.1'. default: none
#' @param ... additional expression directly passed to aes() of ggplot, can refer to colData
#'
#' @return $result: plot, volcano
#'
#' @examples
#' \dontrun{# Volcano plot as overview of results with a result already in 'comp'
#' ... %>%
#' mt_plots_volcano(statname     = "comp",
#'  metab_filter = p.adj < 0.1,
#'  colour       = p.value < 0.05) %>%
#'  ...}
#'
#' @author Jonas Zierer, Jan Krumsiek
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_volcano <- function(D,
                             x = fc,
                             statname,
                             metab_filter = p.value < 0.05,
                             xlabel=gsub("~","",as.character(x)),
                             vline=NA,
                             hline,
                             ggadd=NULL,
                             ...){
    x <- dplyr::enquo(x)

    ## check input
    stopifnot("SummarizedExperiment" %in% class(D))
    if(missing(statname))
        stop("statname must be given for volcanoplot")

    ## rowData
    rd <- rowData(D) %>%
        as.data.frame() %>%
        dplyr::mutate(var = rownames(D))

    ## stat
    data_plot <- mti_get_stat_by_name(D, statname) %>%
        dplyr::inner_join(rd, by = "var") %>%
        dplyr::mutate(xxx = !!x)

    ## SCALE -log10
    reverselog_trans <- function (base = exp(1)){
        trans <- function(x) -log(x, base)
        inv <- function(x) base^(-x)
        scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                          scales::log_breaks(base = base),
                          domain = c(1e-100, Inf))
    }

    ## determine if and where to draw hline
    if (!missing(hline)) {
      hliney <- mti_get_stat_by_name(D, statname) %>%
        dplyr::inner_join(rd, by = "var") %>%
        dplyr::mutate(xxx = !!x) %>% dplyr::filter(!!dplyr::enquo(hline)) %>% .$p.value %>% max()
    } else {
      hliney <- NA
    }

    ## CREATE PLOT
    p <- data_plot %>%
        ## do plot
        ggplot(aes(x = xxx, y = p.value)) +
        # vline?
        (if(!is.na(vline)){geom_vline(xintercept = c(-vline, vline), linetype='dashed', color='#F8766D')}else{NULL}) +
        # hline?
        (if(!is.na(hliney)){geom_hline(yintercept = hliney, linetype='dashed', color='#F8766D')}else{NULL}) +
        # points
        geom_point(aes(...)) +
        scale_y_continuous(trans = reverselog_trans(10),
                           breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        labs(x = xlabel, y = "p-value") +
        ggtitle(statname)

    ## ADD METABOLITE LABELS
    if(!missing(metab_filter)){
        mti_logstatus("add label")
        metab_filter_q <- dplyr::enquo(metab_filter)
        data_annotate <- data_plot %>%
            dplyr::filter(!!metab_filter_q)
        p <- p + ggrepel::geom_text_repel(data = data_annotate,
                           aes(label = name))
    }

    ## ADD AXIS GROUPS
    d <- mti_get_stat_by_name(D, statname, fullstruct=T)
    if ("groups" %in% names(d) && length(d$groups)==2) {
      p <- mti_add_leftright_gg(p, paste0(d$groups[1],' high'), paste0(d$groups[2],' high'))
    }

    # add custom elements?
    if (!is.null(ggadd)) p <- p+ggadd

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
