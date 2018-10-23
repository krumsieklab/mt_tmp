################################################################################
## HEATMAP
################################################################################
#' mt_plot_boxplot
#'
#' create one boxplot per metabolite from SummarizedExperiment and add to metadata
#'
#' @author Jonas Zierer
#' @import SummrizedExperiment
#' @import ggplot2
#' @importFrom dplyr %>% mutate gather left_join filter arrange
#' @importFrom rlang enquo !!
#' @importFrom glue glue
#' @importFrom ggbeeswarm geom_beeswarm
#' @param D SummarizedExperiment object
#' @param stat index of the entry in metadata(D)$results that contains statistic object
#' @param x what phenotype (from colData(D)) should be used on x axis
#' @param ... further parameters forwarded to ggplot::aes
#' @return SummarizedExperiment with boxplots in metadata(D)$results
#' @export mt_plot_boxplot
mt_plots_boxplot <- function(D,
                            stat,
                            x = "x",
                            correct_confounder,
                            metab_filter = p.value < 0.05,
                            metab_sort   = p.value,
                            annotation   = "{sprintf('P-value: %.1e', p.value)}",
                            jitter       = T,
                            rows,
                            cols,
                            ...){
    x <- enquo(x)
    
    if(missing(D) | ! D %of% "SummarizedExperiment")
       stop("D must be a SummarizedExperiment object")

    ## CONFOUNDER
    if(!missing(correct_confounder)){
        message("correcting for ", correct_confounder)
        D <- mti_correctConfounder(D, correct_confounder)
    }
    
    ## rowData
    rd <- rowData(D) %>%
        as.data.frame() %>%
        mutate(var = rownames(D))
    
    ## stat
    if(!missing(stat)){
        if(length(metadata(D)$results) < stat)
            stop("stat element not in results list")
        if("stats" %!in% metadata(D)$results[[stat]]$fun)
            stop("element ", stat, " not a 'stat' object")
        if("var" %!in% colnames(metadata(D)$results[[stat]]$output))
            stop("stat object must have column 'var'")
        stat <- metadata(D)$result[[stat]]$output %>%
                          inner_join(rd, by = "var")
    }else{
        stat <- rd
    }

    ## FILTER METABOLITES
    if(!missing(metab_filter)){
        metab_filter_q <- enquo(metab_filter)
        stat <- stat %>%
            filter(!!metab_filter_q)
        message("filter metabolites: ", metab_filter_q, " [", nrow(stat), " remaining]")
    }
    
    ## SORT METABOLITES
    if(!missing(metab_sort)){
        metab_sort_q <- enquo(metab_sort)
        stat <- stat %>%
            arrange(!!metab_sort_q)
        message("sorted metabolites: ", metab_sort_q)
    }

    ## CREATE PLOT
    p <- D %>%
        mti_format_se_samplewise() %>%
        gather(var, value, one_of(rownames(D))) %>%
        ## add metabolite names
        inner_join(stat, by = "var") %>% 
        ## sort according to stat
        mutate(var = factor(var, levels = stat$var)) %>%
        arrange(var) %>%
        mutate(name = factor(name, levels = unique(name))) %>%
        ## do plot
        ggplot() +
        geom_boxplot(aes(x = !!x, y = value, ...), outlier.shape = ifelse(jitter, NA, 19)) +
        labs(x = NULL, y = NULL)

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
                           x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05 )
    }

    ## SPLIT TO MULTIPLE PAGES
    if(!missing(cols) && !missing(rows)){
        p_plots   <- length(unique(stat$name))
        p_perpage <- cols*rows
        pages     <- ceiling(p_plots / p_perpage)
        ## CREATE PAGES
        message("split ", p_plots, " plots to ", pages, " pages with ", rows, " rows and ", cols, "  cols")
        p <- map(1:pages, ~ p + ggforce::facet_wrap_paginate(~name, scales = "free_y", nrow = rows, ncol  = cols, page = .x))
    }else{
        p <- list(p + facet_wrap(~name, scales = "free_y"))
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






            
