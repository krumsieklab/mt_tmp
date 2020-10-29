#' Plot Box Plots or Scatter Plots
#'
#' Creates one box plot or  scatter plot per metabolite based on given sample annotations
#'
#' @param D \code{SummarizedExperiment} input
#' @param x what phenotype (from colData(D)) should be used on x axis, default "x"
#' @param plot_type "box" or "scatter"
#' @param stat_name index of the entry in metadata(D)$results that contains statistic object
#' @param correct_confounder confounders to adjust for before plotting, formula notation
#' @param metab_filter if given, filter will be applied to data and remaining variables will be labelled in plot, default p.value<0.05
#' @param metab_sort if given, arrange will be applied to data variables will be sorted, default p.value
#' @param annotation if given adds annotation to plot, default = "{sprintf('P-value: %.1e', p.value)}"
#' @param cols number columns of boxplots in $result
#' @param full_info add full information of all sample annotations and statistics results to plottable data.frame? makes plotting more flexible but can render SE objects huge. default: F
#' @param text_size text size of the annotations
#' @param jitter add geom_jitter ("jitter") or geom_beeswarm ("beeswarm") to boxplot, exclude if NULL;  default "beeswarm"
#' @param restrict_to_used_samples whether to filter to the samples that were used in the statistical test, default: T
#' @param manual_ylab manual ylabel (default: none)
#' @param fitline add fit line? (default: T)
#' @param fitline_se add standard error range? (default: T)
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... additional expression directly passed to aes() of ggplot, can refer to colData
#'
#' @return if plot_type = "box", $result: plot, box plot
#' @return if plot_type = "scatter", $result: plot, scatter plot
#'
#' @author JZ, KC
#'
#' @import ggplot2
#'
#' @export


mt_plots_boxplot_scatter <- function(D,
                                     x = "x",
                                     stat_name,
                                     plot_type,
                                     correct_confounder,
                                     metab_filter = p.value < 0.05,
                                     metab_sort = p.value,
                                     annotation = "{sprintf('P-value: %.1e', p.value)}",
                                     cols,
                                     full_info = F,
                                     text_size = 3.88,
                                     jitter = "beeswarm",
                                     restrict_to_used_samples = T,
                                     manual_ylab=NULL,
                                     fitline = T,
                                     fitline_se = T,
                                     ggadd = NULL,
                                     ...){

  stopifnot("SummarizedExperiment" %in% class(D))
  x <- dplyr::enquo(x)

  # create dummy SE so original not changed
  Ds <- D

  ## CONFOUNDER
  if(!missing(correct_confounder)){
    MetaboTools:::mti_logstatus(glue::glue("correcting for {correct_confounder}"))
    Ds <- MetaboTools:::mti_correctConfounder(Ds, correct_confounder)
  }

  ## rowData
  rd <- rowData(Ds) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(Ds))

  ## stat
  if(!missing(stat_name)){
    stat <- MetaboTools:::mti_get_stat_by_name(Ds, stat_name) %>%
      dplyr::inner_join(rd, by = "var")
  }else{
    stat <- rd
    ### KC: ONLY IN BOXPLOT (Why?)
    restrict_to_used_samples <- F # not dependend on a stat
  }

  ## FILTER METABOLITES
  ### KC: metab_filter is never not missing, should there be no default?
  if(!missing(metab_filter)){
    metab_filter_q <- dplyr::enquo(metab_filter)
    stat <- stat %>%
      dplyr::filter(!!metab_filter_q)
    MetaboTools:::mti_logstatus(glue::glue("filter metabolites: {metab_filter_q} [{nrow(stat)} remaining]"))
  }

  ## SORT METABOLITES
  ### KC: metab_sort is never not missing, should there be no default?
  if(!missing(metab_sort)){
    metab_sort_q <- dplyr::enquo(metab_sort)
    stat <- stat %>%
      dplyr::arrange(!!metab_sort_q) %>%
      ## sort according to stat
      dplyr::mutate(name = factor(name, levels = unique(name)))
    MetaboTools:::mti_logstatus(glue::glue("sorted metabolites: {metab_sort_q}"))
  }

  ## CREATE PLOT
  dummy <- Ds %>%
    MetaboTools:::mti_format_se_samplewise() %>% # NOTE: No explosion of dataset size due to active restriction - 6/2/20, JK
    tidyr::gather(var, value, dplyr::one_of(rownames(Ds)))
  ## filter to groups?
  if(plot_type=="box"){
    if (restrict_to_used_samples) {
      filterto <- MetaboTools:::mti_get_stat_by_name(Ds, stat_name, fullstruct=T)$samples.used
      dummy <- dummy[filterto,]
    }
  }

  if(!full_info){
    # filter down only to the variables needed for plotting
    # need to parse x and ... list
    vars <- x %>% dplyr::quo_name()
    q <- dplyr::quos(...)
    if (length(q) > 0) {
      vars <- c(vars, q %>% lapply(function(x){x %>% as.character() %>% gsub("~","",.)}) %>% unlist() %>% as.vector())
    }
    vars <- unique(vars)

    plottitle <- ifelse(missing(stat_name),"",stat_name)
    if(plot_type=="box"){
      # make sure the main outcome variable x is a factor
      mainvar <-x %>% dplyr::quo_name()
      dummy[[mainvar]] <- as.factor(dummy[[mainvar]])

      p <- dummy %>%
        dplyr::select(dplyr::one_of(c("var","value", vars))) %>%
        ## add metabolite names, but only restricted subset from statistics table
        dplyr::inner_join(stat[,dplyr::intersect(colnames(stat),c('var','statistic','p.value','p.adj','name'))], by = "var") %>%
        dplyr::select(-var) %>%
        ## do plot
        ggplot() +
        geom_boxplot(aes(x = as.factor(!!x), y = value, ...), outlier.shape = ifelse(jitter, NA, 19)) +
        labs(x = NULL, y = NULL) +
        ggtitle(plottitle)
    }else{

      p <- dummy %>%
        dplyr::select(dplyr::one_of(c("var","value", vars))) %>%
        ## add metabolite names, but only restricted subset from statistics table
        dplyr::inner_join(stat[,c('var','statistic','p.value','p.adj','name')], by = "var") %>%
        dplyr::select(-var) %>%
        ## do plot
        ggplot() +
        ## add fit line?
        {if (fitline) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fitline_se, color = "black") else NULL} +
        geom_point(aes(x = !!x, y = value, ...)) +
        labs(x = dplyr::quo_name(x), y="metabolite") +
        ggtitle(plottitle)

    }

  }else{

    if(plot_type=="box"){
      # leave full info in
      # can create huge data.frames

      plottitle <- ifelse(missing(stat_name),"",stat_name)
      p <- dummy %>%
        ## add metabolite names
        dplyr::inner_join(stat, by = "var") %>%
        ## do plot
        ggplot() +
        geom_boxplot(aes(x = !!x, y = value, ...), outlier.shape = ifelse(jitter, NA, 19)) +
        labs(x = NULL, y = NULL) +
        ggtitle(plottitle)

    }else{
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
        {if (fitline) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fitline_se, color = "black") else NULL} +
        geom_point(aes(x = !!x, y = value, ...)) +
        labs(x = dplyr::quo_name(x), y="metabolite") +
        ggtitle(plottitle)

    }

  }

  ### BOX PLOT SPECIFIC
  if(plot_type=="box"){
    ## add ylabel
    if (!is.null(manual_ylab)) {
      p <- p + ylab(manual_ylab)
    } else {
      # add label if this is logged data
      r <- Ds %>% MetaboTools:::mti_res_get_path(c("pre","trans","log"))
      if (length(r)>0) {
        p <- p + ylab(r[[1]]$logtxt) # log text contains e.g. "log2"
      }
    }

    ## ADD JITTER
    if(jitter=="beeswarm"){
      p <- p +
        ggbeeswarm::geom_beeswarm(aes(x = !!x, y = value, ...))
    }else if(jitter=="jitter"){
      p <- p +
        geom_jitter(aes(x = !!x, y = value, ...))
    }
  }


  ### COMMON TO BOTH PLOTS
  ## ADD ANNOTATION
  if(!missing(annotation)){
    data_annotate <- stat %>%
      dplyr::mutate(annotate = glue::glue(annotation)) %>%
      dplyr::distinct(name, annotate)
    p <- p + geom_text(data = data_annotate,
                       aes(label = annotate),
                       x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size=text_size )
  }

  if (!is.null(ggadd)) p <- p+ggadd

  ## SPLIT TO MULTIPLE PAGES
  # if there is no plot, create a single empty page
  if (length(unique(stat$name))==0) {
    p <- list(ggplot() + geom_text(aes(x=0, y=0, label='no plots'), size=10))
    output2 <- NULL
  } else {
    p <- p + facet_wrap(.~name, scales = "free_y", ncol=2)
    # fix ggplot environment
    if (D %>% MetaboTools:::mti_get_setting("ggplot_fix")) p <- MetaboTools:::mti_fix_ggplot_env(p)
    p <- list(p)
    output2 <- ceiling(length(unique(stat$name))/2)
  }

  ## add status information & plot
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Metabolite ",ifelse(plot_type=="box", "boxplots", "scatter plots"),", aes: %s", MetaboTools:::mti_dots_to_str(...)),
      output = p,
      output2 = output2
    )
  ## return
  D

}

