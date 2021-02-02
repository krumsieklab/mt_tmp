#' Visualize missing value structure of dataset
#'
#' Creates two kinds of plots:
#' \enumerate{
#'   \item{missing value matrix view (samples and metabolites)}
#'   \item{sorted missing value % plot (samples and/or metabolites)}
#' }
#'
#' @param D \code{SummarizedExperiment} input.
#' @param met_max Which metabolite missingness to mark on the y axis. Default: NA (no line).
#' @param samp_max Which sample missingness to mark on the y axis. Default: NA (no line).
#' @param plot_mets Show metabolite missingness plot? Default: T.
#' @param plot_samples Show sample missingness plot? Default: F.
#' @param sec_axis_mets Add second axis to metabolite missingness plot? Default: F.
#' @param sec_axis_samples Add second axis to sample missingness plot? Default: F.
#' @param sample_labels Which column from colData to use as sample labels? Default: NA.
#' @param plot_data Show entire data missingness plot? Default: T.
#'
#' @return $result$output: plot, two plots
#'
#' @examples
#' \dontrun{%>%  mt_plots_missingness() %>% # without horizontal line
#' %>%  mt_plots_missingness(met_max=0.5) %>% # with horizontal line at 50%
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_missingness <- function(D,
                                 met_max=NA,
                                 samp_max=NA,
                                 plot_mets = T,
                                 plot_samples = F,
                                 sec_axis_mets = F,
                                 sec_axis_samples = F,
                                 sample_labels=NA,
                                 plot_data = T) {

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


