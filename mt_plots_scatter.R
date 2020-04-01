library(SummarizedExperiment)
library(ggplot2)
library(glue)


mt_plots_scatter <- function(D,
                             x = "x",
                             statname,
                             correct_confounder,
                             metab_filter = p.value < 0.05,
                             metab_sort   = p.value,
                             annotation   = "{sprintf('P-value: %.1e', p.value)}",
                             rows,
                             cols,
                             fitline = T,
                             fitline_se = T,
                             full.info=F,
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
  
  
  if (!full.info) { 
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
      ## add fit line?
      {if (fitline) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fitline_se, color = "black") else NULL} + 
      geom_point(aes(x = !!x, y = value, ...)) +
      labs(x = quo_name(x), y="metabolite") +
      ggtitle(plottitle)
    
  } else {
    # leave full info in
    # can create huge data.frames
    
    #
    plottitle <- ifelse(missing(statname),"",statname)
    p <- dummy %>%
      ## add metabolite names
      inner_join(stat, by = "var") %>% 
      ## do plot
      ggplot() +
      ## add fit line?
      {if (fitline) geom_smooth(aes(x = !!x, y = value), method = "lm", se=fitline_se, color = "black") else NULL} + 
      geom_point(aes(x = !!x, y = value, ...)) +
      labs(x = quo_name(x), y="metabolite") +
      ggtitle(plottitle)
  }
  
  
  ## ADD ANNOTATION
  if(!missing(annotation)){
    data_annotate <- stat %>%
      mutate(annotate = glue::glue(annotation)) %>%
      distinct(name, annotate) 
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
      logtxt = sprintf("Metabolite scatter plots, aes: %s", mti_dots_to_str(...)),
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







