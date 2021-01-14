#' Visualize missing value structure of dataset.
#'
#' Creates two kinds of plots:
#' 1. missing value matrix view (samples and metabolites)
#' 2. sorted missing value % plot (samples and/or metabolites)
#'
#' @param D \code{SummarizedExperiment} input
#' @param met_max Which metabolite missingness to mark on the y axis (default: NA = no line)
#' @param samp_max Which sample missingness to mark on the y axis (default: NA = no line)
#' @param plot_mets show metabolite missingness plot? (default: T)
#' @param plot_samples show sample missingness plot? (default: F)
#' @param sec_axis_mets add second axis to metabolite missingness plot? (default: F)
#' @param sec_axis_samples add second axis to sample missingness plot? (default: F)
#' @param sample_labels which column from colData to use as sample labels? (default: NA)
#' @param plot_data show entire data missingness plot? (default: T)
#'
#' @return $result: plot, two plots
#'
#' @examples
#' \dontrun{%>%  mt_plots_missingness() %>% # without horizontal line
#' %>%  mt_plots_missingness(met_max=0.5) %>% # with horizontal line at 50%
#' }
#'
#' @author JK
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_missingness <- function(
  D,         # SummarizedExperiment input
  met_max=NA, # which % to mark on the y axis,
  samp_max=NA, # which % to mark on the y axis,
  plot_mets = T,    # show metabolite missingness plot?
  plot_samples = F, # show sample missingness plot?
  sec_axis_mets = F,
  sec_axis_samples = F,
  sample_labels=NA, # which column from colData to use as sample labels?
  plot_data = T     # show entire data missingess plot?
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
  if (plot_mets) {

    p <-
      data.frame(met=rowData(D)$name, miss=missingness(X)) %>%
      dplyr::arrange(miss) %>%
       ggplot(aes(x=order(miss), met=met)) +
      # geom_point(aes(x=1:ncol(X),y=sort(missingness(X)))) +
      geom_point(aes(y=miss)) +
      xlab("metabolites (sorted)") +
      ylab("missingness") +
      ggtitle("Missing values, metabolites")
    # second axis?
    if(sec_axis_mets){
      p <- p + scale_y_continuous(limits=c(0,1), sec.axis = sec_axis(~.*dim(X)[1], name="# samples")) +
        theme(axis.title.y.right = element_text(angle=90))
    }else{
      p <- p + scale_y_continuous(limits=c(0,1))
    }
    # mark?
    if (!is.na(met_max)){
      p <- p+geom_hline(yintercept=met_max, color='red', linetype="dashed")
    }
    # fix ggplot environment
    if (D %>% mti_get_setting("ggplot_fix")) p <- mti_fix_ggplot_env(p)
    # append
    plots[[length(plots)+1]] <- p
  }

  # missingness %plot per sample
  if (plot_samples) {
    # construct data frame, sample label does not have a default value in MT
    df <- data.frame(miss=missingness(t(X)))
    # add labels?
    if (!is.na(sample_labels)) {
      df$sample = colData(D)[[sample_labels]]
    } else {
      df$sample = rep("", ncol(D))
    }
    # plot
    p <-
      df %>%
      dplyr::arrange(miss) %>%
      ggplot(aes(x=order(miss), sample=sample)) +
      geom_point(aes(y=miss)) +
      #geom_point(aes(x=1:nrow(X),y=sort(missingness(t(X))))) +
      xlab("samples (sorted)") +
      ylab("missingness") +
      ggtitle("Missing values, samples")
    # second axis?
    if(sec_axis_samples){
      p <- p + scale_y_continuous(limits=c(0,1), sec.axis = sec_axis(~.*dim(X)[2], name="# metabolites")) +
        theme(axis.title.y.right = element_text(angle=90))
    }else{
      p <- p + scale_y_continuous(limits=c(0,1))
    }
    # mark?
    if (!is.na(samp_max)){
      p <- p+geom_hline(yintercept=samp_max, color='red', linetype="dashed")
    }
    # fix ggplot environment
    if (D %>% mti_get_setting("ggplot_fix")) p <- mti_fix_ggplot_env(p)
    # append
    plots[[length(plots)+1]] <- p
  }


  if (plot_data) {
    # heatmap
    molten <- reshape2::melt(t(is.na(X)))
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
    if (D %>% mti_get_setting("ggplot_fix")) p <- mti_fix_ggplot_env(p)
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


