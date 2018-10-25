################################################################################
## PCA
################################################################################
#' Principal Component Analysis
#'
#' calculate PCAs and format output
#'
#' @author Jonas Zierer
#' @export pca
#' @param D SummarizedExperiment
#' @param ... additional parameters for prcomp
#' @return list with 1) data.frame containing PCs and other columns from original data
#' @return           2) loadings matrix 
#' @return           3) data.frame contiaining the variance explained by each PC
mt_stats_mutivar_pca <- function(D, fail.on.na=FALSE, ...){
    ass <- assay(D)

    ## EXCLUDE MISSING
    if(any(!is.finite(ass))){
        if(fail.on.na){
            stop("PCA - incomplete samples: exiting")
        }
        nas <- apply(ass, MARGIN = 1, function(x)all(is.finite(x)))
        logmsg(glue::glue("PCA - incomplete metabolites: excluding {sum(!nas)} metabolites - using {sum(nas)}"))
        ass <- ass[nas, ]
    }

    ## EXCLUDE INVARIANT
    var <- apply(ass, MARGIN = 1, function(x)var(x) > 0)
    if(sum(!var) > 0){
        logmsg(glue::glue("PCA - invariant metabolites: excluding {sum(!var)} metabolites - using {sum(var)}"))
        ass <- ass[var, ]
    }

    ## CALCULATE PCA
    pc <- prcomp(t(ass), ...)

    ## VARIANCE EXPLAINED
    lambda.all   <- pc$sdev^2
    data.var.expl <- tibble(var      = colnames(pc$x),
                            num      = 1:length(lambda.all),
                            lambda   = lambda.all,
                            var.expl = lambda.all/sum(lambda.all)) %>%
        mutate(var.expl.cum = cumsum(var.expl)) %>%
        as.data.frame() %>%
        column_to_rownames("var")

    ## LOADINGS
    data.loadings = pc$rotation

    ## X
    data.output <- cbind(as.data.frame(colData(D)),
                         pc$x) 

    ## RETURN
    return(list(pcs = data.output, loadings = data.loadings, explained.variance = data.var.expl))
}


################################################################################
## PCA
################################################################################
#' Plot Principal Component Analysis
#'
#' plot PCAs
#'
#' @author Jonas Zierer
#' @import ggplot2
#' @import SummarizedExperiment
#' @export pca.plot
#' @param D Summarized Experiment
#' @param ... additional aestaetics
#' @return ggplot object
mt_plots_PCA_mult <- function(D,
                         pcs = 3,
                         ol_cut  = NA,
                         ol_col  = NULL,
                         ol_size = 2,
                         ...){

    if(pcs < 2)
        stop("Can't plot less than one PC")
    if(!"SummarizedExperiment" %in% class(D))
        stop("D has to be a SummarizedExperiment")
      
    ## DOPCA
    pca <- mt_stats_mutivar_pca(D)
    
    ## prepare data
    columns    <- str_c("PC", 1:pcs)
    data.pairs <- plotting.makePairs(pca$pcs, cols = str_c("PC", 1:pcs))
    data.plot  <- data.pairs$all %>%
        ## randomize order .. helps to see groups if samples are ordered by phenotye
        arrange(sample(1:n()))
    ## DENSITY
    data.dens  <- data.pairs$densities
    ##
    pc_labeller <- as_labeller(setNames(sprintf("%s (%.1f%%)", columns, pca$explained.variance[columns, "var.expl"]*100),
                                        columns))

    ## PLOT SIMPLE IF ONLY 2 PCS
    if(pcs == 2){
        data.plot <- data.plot %>%
            filter(xvar == "PC1" & yvar == "PC2")
    }
    
    ## PLOT
    p <- ggplot(data.plot, aes(x=x, y=y)) +
        geom_point(aes(...)) +
        facet_grid(yvar ~ xvar, scales = "free", labeller = pc_labeller) +
        ## LABS
        labs(x = NULL, y = NULL) 

    ## ADD DENSITY
    if(pcs > 2){
        p <- p +
        stat_density(aes(x = x, y = ..scaled.. * diff(range(x)) + min(x)), 
                     data = data.dens, position = "identity", 
                     colour = "grey20", geom = "line")
    }

    ## ## LABEL OUTLIERS
    ## if(!is.na(ol_cut) & !is.null(ol_col)){
    ##     data.label <- data.plot %>%
    ##         group_by(xvar, yvar) %>%
    ##         mutate(x_z = (x - mean(x)) / sd(x),
    ##                y_z = (y - mean(y)) / sd(y)) %>%
    ##         filter(abs(x_z) > ol_cut| abs(y_z) > ol_cut) %>%
    ##         ungroup() %>%
    ##         filter(as.character(xvar) > as.character(yvar))
    ##     p.pca <- p.pca +
    ##         geom_text(data = data.label,
    ##                   aes_string(label = ol_col), size = ol_size)
    ## }

    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("PCA, aes: %s", mti_dots_to_str(...)),
                      output = list(p)
                  )
    ## return
    D
}

##############################################################################################################
## MAKE PAIRS FOR SCATTERPLOT MATRIX
##############################################################################################################
## from http://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
#' plotting.makePairs
#'
#' prepare data for scatterplot matrix
#'
#' @author Jonas Zierer
#' @export plotting.makePairs
#' @importFrom dplyr %>% select one_of filter
#' @importFrom tidyr crossing
#' @importFrom purrr map2
#' @param data data containing variables in rows
#' @param cols columns to include in the plot matrix
#' @return dataframe containing each pair of variblaes (former columns)
plotting.makePairs <- function(data, cols = colnames(data)){

    ## SELECT COLUMNS
    cols   <- intersect(cols, colnames(data))
    others <- colnames(data) %without% cols
    
    ## PREPARE GRID
    all <- crossing(xvar = cols, yvar = cols) %>%
        filter(xvar != yvar) %>%
        mutate(data = map2(xvar, yvar, function(xx, yy){
            data %>%
                mutate_(x = xx,
                        y = yy) %>%
                select(x, y, one_of(others))
        })) %>%
        unnest(data) %>%
        mutate(xvar = factor(xvar, levels = cols),
               yvar = factor(yvar, levels = cols))

    ## DENSITIES
    densities <- cols %>%
        map_dfr(function(xx){
            data %>%
               mutate(xvar = xx,
                      yvar = xx) %>%
                mutate_(x = xx) %>%
                select(xvar, yvar, x, one_of(others))
        })

    ## RETURN
    return(list(all=all, densities=densities))
}
