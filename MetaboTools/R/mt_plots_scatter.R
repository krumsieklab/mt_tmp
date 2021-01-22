#' Scatter plots
#'
#' Creates one scatter plot per metabolite based on given sample annotations
#'
#' @param D \code{SummarizedExperiment} input
#' @param x what phenotype (from colData(D)) should be used on x axis, default "x"
#' @param stat_name index of the entry in metadata(D)$results that contains statistic object
#' @param correct_confounder confounders to adjust for before plotting, formula notation
#' @param metab_filter if given, filter will be applied to data and remaining variables will be labelled in plot, default p.value<0.05
#' @param metab_sort if given, arrange will be applied to data variables will be sorted, default p.value
#' @param annotation if given adds annotation to plot, default = "{sprintf('P-value: %.1e', p.value)}"
#' @param rows number rows of boxplots in $result
#' @param cols number columns of boxplots in $result
#' @param fit_line add fit line? (default: T)
#' @param fit_line_se add standard error range? (default: T)
#' @param full_info add full information of all sample annotations and statistics results to plottable data.frame? makes plotting more flexible but can render SE objects huge. default: F
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... additional ggplot parameters
#'
#' @return  $result: plot, scatter plot
#'
#' @examples
#' \dontrun{# scatter plots as overview of results with a result already in 'comp'
#' # color by "age" variable in colData
#' mt_plots_scatter(x                  = age,
#'                  stat_name           = "comp",
#'                  correct_confounder = ~BATCH_MOCK,
#'                  metab_filter       = p.value<0.01,
#'                  metab_sort         = p.value,
#'                  annotation         = "{sprintf('P-value: %.1e', p.value)}\nStatistic: {sprintf('%.2f', statistic)}",
#'                  rows               = 2,
#'                  cols               = 2) %>%
#' }
#'
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_scatter <- function(D,
                             x = "x",
                             stat_name,
                             correct_confounder,
                             metab_filter = p.value < 0.05,
                             metab_sort   = p.value,
                             annotation   = "{sprintf('P-value: %.1e', p.value)}",
                             rows,
                             cols,
                             fit_line = T,
                             fit_line_se = T,
                             full_info=F,
                             ggadd        = NULL,
                             ...){


  stopifnot("SummarizedExperiment" %in% class(D))
  x <- dplyr::enquo(x)

  # create dummy SE so original not changed
  Ds <- D

  ## CONFOUNDER
  if(!missing(correct_confounder)){
    mti_logstatus(glue::glue("correcting for {correct_confounder}"))
    Ds <- mti_correctConfounder(Ds, correct_confounder)
  }

  ## rowData
  rd <- rowData(Ds) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(Ds))

  ## stat
  if(!missing(stat_name)){
    stat <- mtm_get_stat_by_name(Ds, stat_name) %>%
      dplyr::inner_join(rd, by = "var")
  }else{
    stat <- rd
  }

  ## FILTER METABOLITES
  if(!missing(metab_filter)){
    metab_filter_q <- dplyr::enquo(metab_filter)
    stat <- stat %>%
      dplyr::filter(!!metab_filter_q)
    mti_logstatus(glue::glue("filter metabolites: {metab_filter_q} [{nrow(stat)} remaining]"))
  }

  ## SORT METABOLITES
  if(!missing(metab_sort)){
    metab_sort_q <- dplyr::enquo(metab_sort)
    stat <- stat %>%
      dplyr::arrange(!!metab_sort_q) %>%
      ## sort according to stat
      dplyr::mutate(name = factor(name, levels = unique(name)))
    mti_logstatus(glue::glue("sorted metabolites: {metab_sort_q}"))
  }



  ## CREATE PLOT
  dummy <- Ds %>%
    mti_format_se_samplewise() %>% # NOTE: No explosion of dataset size due to active restriction - 6/2/20, JK
    tidyr::gather(var, value, dplyr::one_of(rownames(Ds)))


  if (!full_info) {
    # filter down only to the variables needed for plotting
    # need to parse x and ... list
    vars <- x %>% as.character() %>% gsub("~","",.)
    q <- dplyr::quos(...)
    if (length(q) > 0) {
      vars <- c(vars, q %>% lapply(function(x){x %>% as.character() %>% gsub("~","",.)}) %>% unlist() %>% as.vector())
    }
    vars <- unique(vars)

    #
    plottitle <- ifelse(missing(stat_name),"",stat_name)
    p <- dummy %>%
      dplyr::select(dplyr::one_of(c("var","value", vars))) %>%
      ## add metabolite names, but only restricted subset from statistics table
      dplyr::inner_join(stat[,c('var','statistic','p.value','p.adj','name')], by = "var") %>%
      dplyr::select(-var) %>%
      ## do plot
      ggplot() +
      ## add fit line?
      {if (fit_line) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fit_line_se, color = "black") else NULL} +
      geom_point(aes(x = !!x, y = value, ...)) +
      labs(x = dplyr::quo_name(x), y="metabolite") +
      ggtitle(plottitle)

  } else {
    # leave full info in
    # can create huge data.frames

    #
    plottitle <- ifelse(missing(stat_name),"",stat_name)
    p <- dummy %>%
      ## add metabolite names
      dplyr::inner_join(stat, by = "var") %>%
      ## do plot
      ggplot() +
      ## add fit line?
      {if (fit_line) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fit_line_se, color = "black") else NULL} +
      geom_point(aes(x = !!x, y = value, ...)) +
      labs(x = dplyr::quo_name(x), y="metabolite") +
      ggtitle(plottitle)
  }

  ## ADD ANNOTATION
  if(!missing(annotation)){
    data_annotate <- stat %>%
      dplyr::mutate(annotate = glue::glue(annotation)) %>%
      dplyr::distinct(name, annotate)
    p <- p + geom_text(data = data_annotate,
                       aes(label = annotate),
                       x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05 )
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
      p <- purrr::map(1:pages, ~ p + ggforce::facet_wrap_paginate(~name, scales = "free_y", nrow = rows, ncol  = cols, page = .x))
      ## ADD EMPTY PLOTS TO FILL PAGE
      fill_page <- (pages*p_perpage) - p_plots
      if(fill_page > 0){
        mti_logstatus(glue::glue("add {fill_page} blanks to fill page"))
        spaces <- purrr::map(1:fill_page, ~rep(" ", .x)) %>%
          purrr::map_chr(str_c, collapse = "")
        p[[ pages ]] <- p[[ pages ]] +
          geom_blank(data = data.frame(name = spaces))
      }
    }else{
      p <- p + facet_wrap(.~name, scales = "free_y")
      p <- list(p)
    }
  }

  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Metabolite scatter plots, aes: %s", mti_dots_to_str(...)),
      output = p
    )
  ## return
  D
}

