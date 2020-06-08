#' Visualize missing value structure of dataset.
#' 
#' Creates two plots:
#' 1. missing value matrix view
#' 2. sorted missing value % plot
#'
#' @param D \code{SummarizedExperiment} input 
#' @param metMax Which missingness to mark on the y axis (default: NA = no line)
#'
#' @return $result: plot, two plots
#'
#' @examples
#' %>%  mt_plots_qc_missingness() %>% # without horizontal line
#' %>%  mt_plots_qc_missingness(metMax=0.5) %>% # with horizontal line at 50%
#' 
#' @author JK
#' 
#' @export
mt_plots_qc_missingness <- function(
  D,         # SummarizedExperiment input
  metMax=NA, # which % to mark on the y axis,
  sampMax=NA, # which % to mark on the y axis,
  plot.mets = T,    # show metabolite missingness plot?
  plot.samples = F, # show sample missingness plot?
  sample.labels=NA, # which column from colData to use as sample labels?
  plot.data = T     # show entire data missingess plot?
) {
  
  # helper function
  missingness <- function(X)apply(is.na(X),2,sum)/dim(X)[1]
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))  
  
  # get data
  X <- t(assay(D))
  
  # init
  plots = list()
  
  # missingness %plot per metabolite
  if (plot.mets) {
    
    p <- 
      data.frame(met=rowData(D)$name, miss=missingness(X)) %>% 
      arrange(miss) %>%
       ggplot(aes(x=order(miss), y=miss, met=met)) +
      # geom_point(aes(x=1:ncol(X),y=sort(missingness(X)))) +
      geom_point() +
      ylim(0,1) +
      xlab("metabolites (sorted)") +
      ylab("missingness") +
      ggtitle("Missing values, metabolites")
    # mark?
    if (!is.na(metMax))
      p <- p+geom_hline(yintercept=metMax, color='red', linetype="dashed")
    # fix ggplot environment
    p <- mti_fix_ggplot_env(p)
    # append
    plots[[length(plots)+1]] <- p
  }
  
  # missingness %plot per sample
  if (plot.samples) {
    # construct data frame, sample label does not have a default value in MT
    df <- data.frame(miss=missingness(t(X)))
    # add labels?
    if (!is.na(sample.labels)) {
      df$sample = colData(D)[[sample.labels]]
    } else {
      df$sample = rep("", ncol(D))
    }
    # plot
    p <-
      df %>%
      arrange(miss) %>%
      ggplot(aes(x=order(miss), y=miss, sample=sample)) +
      geom_point() +
      #geom_point(aes(x=1:nrow(X),y=sort(missingness(t(X))))) +
      ylim(0,1) +
      xlab("samples (sorted)") +
      ylab("missingness") +
      ggtitle("Missing values, samples")
    # mark?
    if (!is.na(sampMax))
      p <- p1+geom_hline(yintercept=sampMax, color='red', linetype="dashed")
    # fix ggplot environment
    p <- mti_fix_ggplot_env(p)
    # append
    plots[[length(plots)+1]] <- p
  }
  
  
  if (plot.data) { 
    # heatmap
    molten <- melt(t(is.na(X)))
    colnames(molten) <- c("X1","X2","value") #rename to avoid colname errors 
    molten$X1 <- fixorder(molten$X1)
    molten$X2 <- fixorder(molten$X2)
    p <-
      ggplot(data=molten, aes(x=X1,y=X2)) + 
      geom_tile(aes(fill=value)) +
      theme(axis.text=element_blank(),axis.ticks=element_blank()) + 
      xlab(sprintf("metabolites (%d)", ncol(X))) + ylab(sprintf("samples (%d)", nrow(X))) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      ggtitle("Missing values") + 
      scale_fill_manual(values=c("#FFFFFF", "#000000"), 
                        name="Missing?",
                        breaks=c("FALSE", "TRUE"),
                        labels=c("No", "Yes"))
    # fix ggplot environment
    p <- mti_fix_ggplot_env(p)
    # append
    plots[[length(plots)+1]] <- p
  }
  
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("missingness plots, missing values: %d out of %d (%.2f%%)", sum(is.na(X)), nrow(X)*ncol(X), sum(is.na(X))/(nrow(X)*ncol(X))*100),
      output = plots
    )
  
  # return
  D
  
  
}


