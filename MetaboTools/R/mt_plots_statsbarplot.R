#' mt_plots_statsbarplot
#'
#' Creates a bar plot
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name optional name of the statistics object to be used for filtering
#' @param metab_filter if given, filter will be applied to data and remaining variables will be used to create plot
#' @param aggregate rowData variable used to aggregate variables.
#' @param colorby optional rowData variable used to color barplot. Default NULL.
#' @param yscale plot percentage or frequency of variables. Values c("fraction","count"). Default "percentage".
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... additional expression directly passed to aes() of ggplot, can refer to colData
#'
#' @return $result: plot, barplot
#'
#' @examples
#' \dontrun{# Volcano plot as overview of results with a result already in 'comp'
#' ... %>%
#' mt_plots_statsbarplot(stat_name     = "comp",
#'  metab_filter = p.adj < 0.05,
#'  aggregate    = "SUB_PATHWAY",
#'  colorby      = "SUPER_PATHWAY",
#'  yscale       = "count") %>%
#'  ...}
#'
#' @author Elisa Benedetti
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_statsbarplot <- function(D,
                             stat_name,
                             metab_filter = p.value < 0.05,
                             aggregate = "SUB_PATHWAY",
                             colorby = NULL,
                             ggadd = NULL,
                             yscale = "fraction",
                             ...){

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_name) & !missing(metab_filter))
    stop("stat_name must be given for metab_filter to work.")
  if(!is.null(colorby))
    if(!(colorby %in% colnames(rowData(D))))
    stop("colour is not in rowData")

  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))

  perc <- rd[[aggregate]] %>%
    unlist %>% table %>% as.data.frame()
  colnames(perc) <- c("name","count")
  ## subselect variables
  if(!missing(metab_filter)) {
    metab_filter_q <- dplyr::enquo(metab_filter)
    sel <- MetaboTools:::mti_get_stat_by_name(D=D,name=stat_name) %>%
    dplyr::filter(!!metab_filter_q) %>%
      dplyr::select(var)
    rd <- rd %>%
      dplyr::filter(var %in% sel$var)
  }
  # throw an error if filtering gives an empty matrix
  if(nrow(rd)==0)
    stop("The filtering step gave an empty matrix.")
  
  ## get data to plot
  data_plot <- rd[[aggregate]] %>% 
    unlist %>% table %>% as.data.frame()
  colnames(data_plot) <- c("name","count")
  # add number of metabolites in each pathway
  perc %<>% dplyr::filter(name %in% data_plot$name) 
  perc <- perc[match(data_plot$name,perc$name),]
  data_plot <- data_plot %>%
    dplyr::mutate(label=sprintf("%s [%d]", name, perc$count)) %>%
    # add fraction variable
    dplyr::mutate(fraction= count/perc$count)
    
  # add color column if not given
  if(is.null(colorby)) {
    colorby <- paste(aggregate,"color", collapse = "_")
    rd[[colorby]] <- "pathway"
  } 
  
  # create dictionary between aggregate and colorby variables
  dict <- rd %>% dplyr::select(!!sym(aggregate),!!sym(colorby)) %>% tidyr::unnest(cols=c(aggregate)) %>% as.data.frame()
  dict <- dict[!duplicated(dict[[aggregate]]),]
  
  # add color to data_plot
  data_plot <- data_plot %>%
    dplyr::left_join(dict, by=c("name"=aggregate)) %>%
    dplyr::rename(color=sym(colorby))
  # convert labels to factor for sorting
  data_plot$label <- as.factor(data_plot$label)
  
  ## CREATE PLOT
  p <- data_plot %>%
    ggplot(aes(x=label, y=!!sym(yscale), fill=color)) + 
    geom_bar(stat = "identity") +
    (if(yscale=="fraction") {ggtitle(sprintf("Fraction of pathway affected, %s",  gsub("~", "", rlang::expr_text(enquo(metab_filter)))))}else{ggtitle(sprintf("Number of hits per pathway, %s",  gsub("~", "", rlang::expr_text(enquo(metab_filter)))))}) +
    labs(x="",fill = colorby) +
    scale_x_discrete(limits = rev(levels(data_plot$label))) +
    coord_flip()

  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  # fix ggplot environment
  p <- MetaboTools:::mti_fix_ggplot_env(p)

  ## add status information & plot
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = ifelse(exists("stat_name"), sprintf("bar plot for comparison %s, by %s, filtered for %s, using %s", stat_name, aggregate, gsub("~", "", rlang::expr_text(enquo(metab_filter))), yscale),
                      sprintf("bar plot by %s using %s", aggregate, yscale)),
      output = list(p)
    )
  ## return
  D
}

