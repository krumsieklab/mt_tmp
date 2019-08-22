require(SummarizedExperiment)
require(ggplot2)
require(glue)

#' mt_plot_boxplot
#'
#' Creates one boxplot per metabolite from SummarizedExperiment and add to metadata
#'
#' @param D \code{SummarizedExperiment} input
#' @param x what phenotype (from colData(D)) should be used on x axis, default "x"
#' @param statname index of the entry in metadata(D)$results that contains statistic object
#' @param correct_confounder confounders to adjust for before plotting
#' @param metab_filter if given, filter will be applied to data and remaining variables will be labelled in plot, default p.value<0.05
#' @param metab_sort if given, arrange will be applied to data variables will be sorted, default p.value
#' @param annotation if given adds annotation to plot, default = "{sprintf('P-value: %.1e', p.value)}",
#' @param text.size text size of the annotations
#' @param jitter whether to add jitter to boxplot,  default T
#' @param rows number rows of boxplots in $result
#' @param cols number columns of boxplots in $result
#' @param restrict_to_groups whether to filter by groups, default T
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... additional expression directly passed to aes() of ggplot, can refer to colData
#' 
#' @return  $result: plot, boxplot
#' 
#' @examples 
#' # boxplots as overview of results with a result already in 'comp'
#' # color by "Group" variable in colData
#' mt_plots_boxplot(x                  = Group,
#'                  statname           = "comp",
#'                  correct_confounder = ~BATCH_MOCK,
#'                  metab_filter       = p.value<0.01,
#'                  metab_sort         = p.value,
#'                  annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
#'                  rows               = 2,
#'                  cols               = 2) %>%
#'                  ...
#'                  
#'  @author Jonas Zierer
mt_plots_boxplot <- function(D,
                             x = "x",
                             statname,
                             correct_confounder,
                             metab_filter = p.value < 0.05,
                             metab_sort   = p.value,
                             annotation   = "{sprintf('P-value: %.1e', p.value)}",
                             text.size = 3.88,
                             jitter       = T,
                             rows,
                             cols,
                             restrict_to_groups=T,
                             ggadd        = NULL,
                             ...){
  
  stopifnot("SummarizedExperiment" %in% class(D))
  x <- enquo(x)
  
  
  ## CONFOUNDER
  if(!missing(correct_confounder)){
    mti_logstatus(glue::glue("correcting for {correct_confounder}"))
    D <- mti_correctConfounder(D, correct_confounder)
  }
  
  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    mutate(var = rownames(D))
  
  ## stat
  if(!missing(statname)){
    stat <- mti_get_stat_by_name(D, statname) %>%
      inner_join(rd, by = "var")
  }else{
    stat <- rd
    restrict_to_groups <- F # not dependend on a stat
  }
  
  ## FILTER METABOLITES
  if(!missing(metab_filter)){
    metab_filter_q <- enquo(metab_filter)
    stat <- stat %>%
      filter(!!metab_filter_q)
    mti_logstatus(glue::glue("filter metabolites: {metab_filter_q} [{nrow(stat)} remaining]"))
  }
  
  ## SORT METABOLITES
  if(!missing(metab_sort)){
    metab_sort_q <- enquo(metab_sort)
    stat <- stat %>%
      arrange(!!metab_sort_q) %>%
      ## sort according to stat
      mutate(name = factor(name, levels = unique(name)))
    mti_logstatus(glue::glue("sorted metabolites: {metab_sort_q}"))
  }
  
  
  
  ## CREATE PLOT
  dummy <- D %>%
    mti_format_se_samplewise() %>%
    gather(var, value, one_of(rownames(D)))
  ## filter to groups?
  if (restrict_to_groups) {
    filterto <- mti_get_stat_by_name(D, statname, fullstruct=T)$groups
    dummy <- dummy[dummy[[stat$term[1]]] %in% filterto,]
  }
  
  

  # filter down only to the variables needed for plotting
  # need to parse x and ... list
  vars <- x %>% as.character() %>% gsub("~","",.)
  q <- quos(...)
  if (length(q) > 0) {
    vars <- c(vars, q %>% lapply(function(x){x %>% as.character() %>% gsub("~","",.)}) %>% unlist() %>% as.vector())
  }
  vars <- unique(vars)
  
  # make sure the main outcome variable x is a factor
  mainvar <- x %>% as.character() %>% gsub("~","",.)
  dummy[[mainvar]] <- as.factor(dummy[[mainvar]])
  
  #
  plottitle <- ifelse(missing(statname),"",statname)
  p <- dummy %>%
    dplyr::select(one_of(c("var","value", vars))) %>%
    ## add metabolite names, but only restricted subset from statistics table
    inner_join(stat[,c('var','statistic','p.value','p.adj','name')], by = "var") %>% 
    dplyr::select(-var) %>%
    ## do plot
    ggplot() +
    geom_boxplot(aes(x = as.factor(!!x), y = value, ...), outlier.shape = ifelse(jitter, NA, 19)) +
    labs(x = NULL, y = NULL) +
    ggtitle(plottitle)
  
  ## add ylabel if this is logged data
  r <- D %>% mti_res_get_path(c("pre","trans","log"))
  if (length(r)>0) {
    p <- p + ylab(r[[1]]$logtxt) # log text contains e.g. "log2"
  }
  
  ## ADD JITTER
  if(jitter){
    p <- p +
      ggbeeswarm::geom_beeswarm(aes(x = !!x, y = value, ...))
  }
  
  ## ADD ANNOTATION
  if(!missing(annotation)){
    data_annotate <- stat %>%
      mutate(annotate = glue::glue(annotation)) %>%
      distinct(name, annotate) 
    p <- p + geom_text(data = data_annotate,
                       aes(label = annotate),
                       x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size=text.size )
  }
  
  if (!is.null(ggadd)) p <- p+ggadd
  
  ## SPLIT TO MULTIPLE PAGES
  # if there is no plot, create a single empty page
  if (length(unique(stat$name))==0) {
    p <- list(ggplot() + geom_text(aes(x=0, y=0, label='no plots'), size=10))
  } else {
    if(!missing(cols) && !missing(rows)){
      p_plots   <- length(unique(stat$name))
      p_perpage <- cols*rows
      pages     <- ceiling(p_plots / p_perpage)
      mti_logstatus(glue::glue("split {p_plots} plots to {pages} pages with {rows} rows and {cols} cols"))
      p <- map(1:pages, ~ p + ggforce::facet_wrap_paginate(~name, scales = "free_y", nrow = rows, ncol  = cols, page = .x))
      ## ADD EMPTY PLOTS TO FILL PAGE
      fill_page <- (pages*p_perpage) - p_plots
      if(fill_page > 0){
        mti_logstatus(glue::glue("add {fill_page} blanks to fill page"))
        spaces <- map(1:fill_page, ~rep(" ", .x)) %>%
          map_chr(str_c, collapse = "")
        p[[ pages ]] <- p[[ pages ]] +
          geom_blank(data = data.frame(name = spaces))
      }
    }else{
      p <- list(p + facet_wrap(.~name, scales = "free_y"))
    }
  }
  
  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Metabolite boxplot, aes: %s", mti_dots_to_str(...)),
      output = p
    )
  ## return
  D
}


mti_correctConfounder <- function(D, formula){
  d <- D %>% mti_format_se_samplewise()
  d_cor <- rownames(D) %>%
    map_dfc(function(m){
      f   <- update.formula(formula, str_c(m, "~."))
      mod <- lm(f, data = d, na.action = na.exclude)
      res <- resid(mod)
      res
    }) %>%
    setNames(rownames(D)) %>%
    as.matrix() %>% t()
  colnames(d_cor) <- colnames(D)
  assay(D)        <- d_cor
  D
}







