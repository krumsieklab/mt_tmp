# MetaboTools
#
# Visualize missing value structure of dataset.
#
# Creates two plots:
# 1. missing value matrix view
# 2. sorted missing value % plot
#
# last update: 2018-10-24
# authors: JK
#

# helper function
missingness <- function(X)apply(is.na(X),2,sum)/dim(X)[1]

# main function
mt_plots_qc_missingness <- function(
  D # SummarizedExperiment input
) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))  
  
  # get data
  X <- t(assay(D))
  
  # missingness % plot
  p1 <- ggplot() +
    geom_point(aes(x=1:ncol(X),y=sort(missingness(X)))) +
    xlab("metabolites (sorted)") +
    ylab("missingness") +
    ggtitle("Missing values")
  # heatmap
  molten <- melt(t(is.na(X)))
  colnames(molten) <- c("X1","X2","value") #rename to avoid colname errors 
  molten$X1 <- fixorder(molten$X1)
  molten$X2 <- fixorder(molten$X2)
  p2 <-
    ggplot(data=molten, aes(x=X1,y=X2)) + 
    geom_tile(aes(fill=value)) +
    theme(axis.text=element_blank(),axis.ticks=element_blank()) + 
    xlab("metabolites") + ylab("samples") +
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
      logtxt = "missingness plots",
      output = list(p1,p2)
    )
  
  # return
  D
  
  
}