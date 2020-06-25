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
#' @param sort sort pathways in plot according to yscale. Default FALSE
#' @param assoc_sign optional parameter to discriminate between positive and negative associations.
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
#'  yscale       = "count",
#'  sort         = TRUE) %>%
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
                                  sort = FALSE,
                                  assoc_sign,
                                  ...){
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_name) & !missing(metab_filter))
    stop("stat_name must be given for metab_filter to work.")
  if(missing(stat_name) & !missing(assoc_sign))
    stop("stat_name must be given for assoc_sign to work.")
  if(!is.null(colorby))
    if(!(colorby %in% colnames(rowData(D))))
      stop("colour is not in rowData")
  
  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))
  
  perc <- rd[[aggregate]] %>%
    unlist %>% table(exclude = NULL) %>% as.data.frame()
  colnames(perc) <- c("name","count")
  
  ## subselect variables
  if(!missing(metab_filter)) {
    metab_filter_q <- dplyr::enquo(metab_filter)
    sel <- MetaboTools:::mti_get_stat_by_name(D=D,name=stat_name) %>%
      dplyr::filter(!!metab_filter_q)
    rd <- rd %>%
      dplyr::filter(var %in% sel$var)
  }
  
  # if filtering gives an empty matrix, produce an empty plot
  if(nrow(rd)==0) {
    
    p <- ggplot() + 
      geom_text(aes(x=0,y=0, label="No significant results"), size=10)
    
  } else {
    
    if(!missing(assoc_sign)){
      if(!(assoc_sign %in% colnames(sel))) {
        stop(sprintf("Could not find column called %s in the statistical results called %s", assoc_sign, stat_name))
      } else {
        # reorder sel according to rd
        sel <- sel[match(sel$var,rd$var),] %>%
          dplyr::mutate(association=ifelse(sign(!!sym(assoc_sign))>0, "positive", "negative")) %>%
          dplyr::select(var,association)
        
        data_plot <- cbind.data.frame(name=rd[[aggregate]] %>% unlist, 
                                      association=rep(sel$association, times= (rd[[aggregate]] %>% sapply(length)))) %>% 
          table(exclude = NULL) %>% as.data.frame()
        colnames(data_plot) <- c("name","association","count")
        
      } 
    } else {
      sel <- sel[match(sel$var,rd$var),] %>%
        dplyr::select(var)
      
      data_plot <- rd[[aggregate]] %>% 
        unlist %>% table(exclude = NULL) %>% as.data.frame()
      colnames(data_plot) <- c("name","count")
    }
    
    # add number of metabolites in each pathway
    perc <- data_plot %>% dplyr::select(name) %>%
      left_join(perc, by="name")
    data_plot <- data_plot %>%
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
    
    # if pathway mapping exists in the metadata, use the names provided there
    x <- D %>% metadata
    if ("pathways" %in% names(x)){
      if (aggregate %in% names(x$pathways)) {
        # add pathway names to dataframe
        data_plot %<>% dplyr::left_join(x$pathways[[aggregate]][,c("ID","pathway_name")], by=c("name"="ID"))
        # substitute codes for names
        data_plot$name <- data_plot$pathway_name
      } else{
        warning(sprintf("%s field not found in the metadata",aggregate))
      } 
    }
    # create labels for plotting
    data_plot %<>% dplyr::mutate(label=sprintf("%s [%d]", name, perc$count))
    
    # convert labels to factor to sort alphabetically
    data_plot$label <- as.factor(data_plot$label)
    # optional sorting
    if (sort){
      data_plot$label <- reorder(data_plot$label, -data_plot[[yscale]])
    } 
    
    ## CREATE PLOT
    p <- data_plot %>%
      ggplot() + 
      geom_bar(aes(x=label, y=!!sym(yscale), fill=color), stat = "identity") +
      (if("association" %in% colnames(data_plot)) {geom_bar(aes(x=label, y=!!sym(yscale), fill=color, color=association), size=1, stat = "identity")}) +
      # (if("association" %in% colnames(data_plot)) {geom_bar(aes(x=label, y=!!sym(yscale), fill=color, color=association, size=2), stat = "identity")}else{geom_bar(aes(x=label, y=!!sym(yscale), fill=color), stat = "identity")}) +
      scale_color_manual(values=c("green", "blue")) +
      (if(yscale=="fraction") {ggtitle(sprintf("Fraction of pathway affected, %s",  gsub("~", "", rlang::expr_text(enquo(metab_filter)))))}else{ggtitle(sprintf("Number of hits per pathway, %s",  gsub("~", "", rlang::expr_text(enquo(metab_filter)))))}) +
      labs(x="",fill = colorby) +
      scale_x_discrete(limits = rev(levels(data_plot$label))) +
      coord_flip()
    
    # add custom elements?
    if (!is.null(ggadd)) p <- p+ggadd
  }
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

