#' mt_plots_statsbarplot
#'
#' Creates a bar plot
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name list of names of the statistics objects to be used for filtering
#' @param metab_filter filter will be applied to data and remaining variables will be used to create plot
#' @param aggregate rowData variable used to aggregate variables.
#' @param colorby optional rowData variable used to color barplot. Default NULL.
#' @param yscale plot percentage or frequency of variables. Values c("fraction","count"). Default "fraction".
#' @param sort sort pathways in plot according to yscale. Default FALSE.
#' @param assoc_sign optional parameter to discriminate between positive and negative associations. Needs to be the name of a column in the statistical results indicated by stat_name.
#' @param keep.unmapped boolean, if TRUE keeps metabolites with no pathway annotations. Default "FALSE".
#' @param add_empty boolean, if TRUE adds also empty pathways to the barplot.
#' @param output.file optional Excel filename to save data to 
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... additional expression directly passed to aes() of ggplot, can refer to colData
#'
#' @return $result: plot, barplot
#'
#' @examples
#' \dontrun{# Barplot as overview of results with a result already in 'comp'
#' ... %>%
#' mt_plots_statsbarplot(stat_name     = "comp",
#'  metab_filter = p.adj < 0.05,
#'  aggregate    = "SUB_PATHWAY",
#'  colorby      = "SUPER_PATHWAYdevtools::install(codes.makepath("MT/MetaboTools"))",
#'  yscale       = "count",
#'  sort         = TRUE) %>%
#'  ...}
#'
#' @author Elisa Benedetti
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import openxlsx
#' @import SummarizedExperiment
#'
#' @export

mt_plots_statsbarplot <- function(D,
                                  stat_name,
                                  metab_filter = p.value < 1,
                                  aggregate = "SUB_PATHWAY",
                                  colorby = NULL,
                                  ggadd = NULL,
                                  yscale = "fraction",
                                  sort = FALSE,
                                  assoc_sign,
                                  add_empty = FALSE,
                                  keep.unmapped = FALSE,
                                  output.file = NULL,
                                  ...){

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_name) & !missing(metab_filter))
    stop("stat_name must be given for metab_filter to work.")
  if(missing(stat_name) & !missing(assoc_sign))
    stop("stat_name must be given for assoc_sign to work.")
  if(!(aggregate %in% colnames(rowData(D))))
    stop(sprintf("aggregate column '%s' not found in rowData", aggregate))
  if(!is.null(colorby))
    if(!(colorby %in% colnames(rowData(D))))
      stop(sprintf("colorby column '%s' not found in rowData", colorby))

  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))
  # set the nulls to unknown
  if(keep.unmapped){
    rd[[aggregate]][which(rd[[aggregate]]=="NULL")] <- "Unmapped"
  }

  perc <- rd[[aggregate]] %>%
    unlist %>% table(exclude = NULL) %>% as.data.frame()
  colnames(perc) <- c("name","count")

  flag_filter <- ifelse((!missing(metab_filter)), T,F)
  flag_sign <- ifelse((!missing(assoc_sign)), T,F)

  data <- lapply(stat_name %>% {names(.)=.;.}, function(ss){
    ## subselect variables
    if(flag_filter) {
      metab_filter_q <- dplyr::enquo(metab_filter)
      sel <- MetaboTools:::mti_get_stat_by_name(D=D,name=ss) %>%
        dplyr::filter(!!metab_filter_q)
      rd <- rd %>%
        dplyr::filter(var %in% sel$var)

    }

    # if filtering gives an empty matrix, produce an empty df
    if(nrow(rd)==0) {
      data_plot <- data.frame(name=as.character(),
                              count = as.numeric())
      anno <- data.frame()
      
    } else {
      # if assoc_sign given, include in data
      if(flag_sign){
        if(!(assoc_sign %in% colnames(sel))) {
          stop(sprintf("Could not find column called %s in the statistical results called %s", assoc_sign, stat_name))
        } else {
          # reorder sel according to rd
          sel <- sel[match(sel$var,rd$var),] %>%
            dplyr::mutate(association=ifelse(sign(!!sym(assoc_sign))>0, "positive", "negative")) %>%
            dplyr::select(var,association)
          
          # create data.frame for plotting 
          data_plot <- data.frame(name=rd[[aggregate]] %>% unlist %>% as.vector,
                                  association=rep(sel$association, times= (rd[[aggregate]] %>% sapply(length)))) %>%
            table(exclude = NULL) %>% as.data.frame()
          colnames(data_plot) <- c("name","association","count")

        }
      } else {
        # reorder sel according to rd
        sel <- sel[match(sel$var,rd$var),] %>%
          dplyr::select(var)
        
        # create data.frame for plotting
        data_plot <- rd[[aggregate]] %>%
          unlist %>% table(exclude = NULL) %>% as.data.frame()
        colnames(data_plot) <- c("name","count")

      }

      if(add_empty){
        # check which aggregate entries are not included
        agg <- rowData(D) %>% as.data.frame() %>% .[[aggregate]] %>% unlist %>% unique
        agg_empty <- agg[which(!(agg %in% unique(as.character(data_plot$name))))]

        # create data frame
        empty <- data.frame(name = agg_empty, count = rep(0, times=length(agg_empty)))
        
        if("association" %in% colnames(data_plot)){
          empty$association <- "positive"
        }

        # add to data
        data_plot %<>% dplyr::full_join(empty, by=colnames(data_plot))
      }

      # add number of metabolites in each pathway
      perc <- data_plot %>% dplyr::select(name) %>%
        dplyr::left_join(perc, by="name")
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
      
      # create annotation data
      anno <- anno <- data.frame(name = rep(rd$name, times=sapply(rd[[aggregate]], length) %>% as.vector()), 
                                 var = rep(rd$var, times=sapply(rd[[aggregate]], length) %>% as.vector()), 
                                 pathway = unlist(rd[[aggregate]]),
                                 color = ifelse(!is.null(colorby), rep(rd[[colorby]], times=sapply(rd[[aggregate]], length) %>% as.vector()),"pathway")) %>%
        dplyr::left_join(MetaboTools:::mti_get_stat_by_name(D=D,name=ss) , by="var") %>%
        dplyr::select(-var)

      # if pathway mapping exists in the metadata, use the names provided there
      x <- D %>% metadata
      if ("pathways" %in% names(x)){
        if (aggregate %in% names(x$pathways)) {
          # add pathway names to dataframe
          data_plot %<>% 
            dplyr::left_join(x$pathways[[aggregate]][,c("ID","pathway_name")], by=c("name"="ID"))
          anno %<>% 
            dplyr::left_join(x$pathways[[aggregate]][,c("ID","pathway_name")], by=c("pathway"="ID")) %>%
            dplyr::select(name,pathway,pathway_name,color,everything())
          
          # set Unknown pathway names to Unknown
          if(length(which(is.na(data_plot$pathway_name)))>0){
            data_plot$pathway_name[which(is.na(data_plot$pathway_name))] <- "Unknown"
          }
          # substitute codes for names
          data_plot$name <- data_plot$pathway_name
        } else{
          warning(sprintf("%s field not found in the metadata", aggregate))
        }
      }
      # create labels for plotting
      data_plot %<>% dplyr::mutate(label=sprintf("%s [%d]", name, perc$count))

      # convert labels to factor to sort alphabetically
      data_plot$label <- as.factor(data_plot$label)

      # add comparison name to df
      data_plot$comp <- ss
      anno$comp <- ss

    }
    list(dt = data_plot, anno = anno)
  })
  
  # function to revert string structure
  revert_list_str_4 <- function(ls) {
    # get sub-elements in same order
    x <- lapply(ls, `[`, names(ls[[1]]))
    # stack and reslice
    apply(do.call(rbind, x), 2, as.list) 
  }
  
  data <- revert_list_str_4(data)

  # if there is at least one result, produce plot, otherwise output empty plot
  if((sapply(data$dt, function(ss){dim(ss)[1]}) %>% sum()) >0) {
    # merge list into a single dataframe
    data_plot <- do.call(rbind, data$dt) %>% as.data.frame()
    anno <- do.call(rbind, data$anno) %>% as.data.frame()

    # optional sorting (only for single statistical results)
    if (sort){
      data_plot$label <- reorder(data_plot$label, -data_plot[[yscale]])
    }
    # sort comp so that facets appear in the same order given by the user
    data_plot$comp <- factor(data_plot$comp,levels=stat_name)

    # convert count to numeric
    data_plot$count %<>% as.numeric

    ## CREATE PLOT
    p <- ggplot(data_plot, aes(label)) +
      (if("association" %in% colnames(data_plot)) {geom_bar(data = subset(data_plot, association == "positive"), aes(y = !!sym(yscale), fill = color), stat = "identity", position = "dodge", color="black", size=0.4)}) +
      (if("association" %in% colnames(data_plot)) {geom_bar(data = subset(data_plot, association == "negative"), aes(y = -!!sym(yscale), fill = color), stat = "identity", position = "dodge", color="black", size=0.4)} else{geom_bar(aes(x=label, y=!!sym(yscale), fill=color), stat = "identity", color="black", size=0.4)}) +
      (if(yscale=="fraction") {ggtitle(sprintf("Fraction of pathway affected, %s", gsub("~", "", rlang::expr_text(enquo(metab_filter)))))}else{ggtitle(sprintf("Number of hits per pathway, %s", gsub("~", "", rlang::expr_text(enquo(metab_filter)))))}) +
      (if(yscale=="count" & "association" %in% colnames(data_plot)) {expand_limits(y=c(-max(data_plot$count, na.rm = T)*1.7, max(data_plot$count, na.rm = T)*1.7))}) +
      (if(yscale=="count" & !("association" %in% colnames(data_plot))) {expand_limits(y=c(0, max(data_plot$count, na.rm = T)*1.7))}) +
      (if(yscale=="fraction" & "association" %in% colnames(data_plot)) {expand_limits(y=c(-1, 1))}) +
      geom_hline(yintercept = 0,colour = "black", size=0.4) +
      labs(x="",fill = colorby) +
      theme(plot.title = element_text(hjust = 0.4)) +
      scale_x_discrete(limits = rev(levels(data_plot$label)))

    # add phenotype labels to x axis
    if("association" %in% colnames(data_plot) & length(stat_name)==1){
      d <- MetaboTools:::mti_get_stat_by_name(D, stat_name, fullstruct=T)
      if ("groups" %in% names(d) && length(d$groups)==2) {
        # get breaks
        ggbld <- ggplot2::ggplot_build(p)
        yticks = ggbld$layout$panel_params[[1]]$y$minor_breaks # using minor_breaks because sometimes breaks would not work
        # edit labels to include groups
        ytlabs = yticks
        ytlabs[1] <- sprintf("%s\n%s", yticks[1], sprintf("high in %s", d$groups[1]))
        ytlabs[length(ytlabs)] <- sprintf("%s\n%s", yticks[length(yticks)], sprintf("high in %s", d$groups[2]))
        # apply new labels
        p <- p +
          scale_y_continuous(breaks = yticks, labels = ytlabs)
      }
    }

    # flip axes and add annotations on bars
    p <- p +
      coord_flip() +
      (if(yscale=="count" & !("association" %in% colnames(data_plot))) {geom_text(data=data_plot, aes(label, !!sym(yscale), label= sprintf("%.2f%%", fraction*100)),
                                                                               position = position_dodge(width=0.9), hjust = -0.1, size=2.5)}) +
      (if(yscale=="count" & "association" %in% colnames(data_plot)) {geom_text(data=dplyr::filter(data_plot, association=="positive"), aes(label, !!sym(yscale), group= association,label= sprintf("%.2f%%", fraction*100)),
                position = position_dodge(width=0.9), hjust = -0.1, size=2.5)}) +
      (if(yscale=="count" & "association" %in% colnames(data_plot)) {geom_text(data=data_plot %>% dplyr::filter(association=="negative"), aes(label, -!!sym(yscale), group= association,label= sprintf("%.2f%%", fraction*100)),
                position = position_dodge(width=0.9), hjust = 1.1, size=2.5)}) +
      facet_wrap(~comp)

    # add custom elements?
    if (!is.null(ggadd)) p <- p + ggadd
    
    # save plot parameters to be passed to the html generator for dynamical plot height
    re <- p %>%
      ggplot2::ggplot_build() %>%
      magrittr::extract2('layout') %>% 
      magrittr::extract2('layout')

    nr <- data_plot$name %>% unique %>% length # number of pathways
    ncol <- re$COL %>% max() # number of panel columns
    nrow <- re$ROW %>% max() # number of panel rows

    # fix ggplot environment
    if (D %>% MetaboTools:::mti_get_setting("ggplot_fix")) p <- MetaboTools:::mti_fix_ggplot_env(p)

  } else {

    p <- ggplot() +
      geom_text(aes(x=0,y=0, label="No significant results"), size=10)
    
    # save plot parameters to be passed to the html generator for dynamical plot height
    nr <- 0 # number of pathways
    ncol <- NULL
    nrow <- NULL

  }

  if(!is.null(output.file)){
    if(exists("data_plot")){
      wb = createWorkbook()
      sheet = addWorksheet(wb, "Parameters")
      writeData(wb, sheet=sheet, list(comparisons = stat_name, metab_filter = gsub("~", "", rlang::expr_text(enquo(metab_filter))), aggregate = aggregate, coloredby = colorby))
      sheet = addWorksheet(wb, "AggregatedPathways")
      writeData(wb, sheet=sheet, data_plot, rowNames = F, colNames = T)
      sheet = addWorksheet(wb, "IndividualResults")
      writeData(wb, sheet=sheet, anno, rowNames = F, colNames = T)
      saveWorkbook(wb, output.file, overwrite = T)
    } else {
      warning("mt_plots_statsbarplot: No significant results. output.file ignored.")
    }
  }

  ## add status information & plot
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = ifelse(exists("stat_name"), sprintf("bar plot for comparison %s, by %s, filtered for %s, using %s", paste(stat_name,collapse = ", "), aggregate, gsub("~", "", rlang::expr_text(enquo(metab_filter))), yscale),
                      sprintf("bar plot by %s using %s", aggregate, yscale)),
      output = list(p),
      output2 = list(nr = nr, npancol = ncol, npanrow = nrow)
    )
  ## return
  D
}

