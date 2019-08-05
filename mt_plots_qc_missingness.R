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
mt_plots_qc_missingness <- function(
  D,        # SummarizedExperiment input
  metMax=NA # which % to mark on the y axis
) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))  
  
  # get data
  X <- t(assay(D))
  
  # missingness % plot
  p1 <- ggplot() +
    geom_point(aes(x=1:ncol(X),y=sort(missingness(X)))) +
    ylim(0,1) +
    xlab("metabolites (sorted)") +
    ylab("missingness") +
    ggtitle("Missing values")
  # mark?
  if (!is.na(metMax))
    p1 <- p1+geom_hline(yintercept=metMax, color='red', linetype="dashed")
  
  # heatmap
  molten <- melt(t(is.na(X)))
  colnames(molten) <- c("X1","X2","value") #rename to avoid colname errors 
  molten$X1 <- fixorder(molten$X1)
  molten$X2 <- fixorder(molten$X2)
  p2 <-
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
  

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("missingness plots, missing values: %d out of %d (%.2f%%)", sum(is.na(X)), nrow(X)*ncol(X), sum(is.na(X))/(nrow(X)*ncol(X))*100),
      output = list(p1,p2)
    )
  
  # return
  D
  
  
}

# helper function
missingness <- function(X)apply(is.na(X),2,sum)/dim(X)[1]
