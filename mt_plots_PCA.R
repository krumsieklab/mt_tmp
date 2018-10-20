# MetaboTools
#
# PCA plot, show first two PCs plotted against each other.
# - Coloring by factor or numerical.
# - Symbols by factor.
#
# last update: 2018-10-13
# authors: JK
#


require(ggplot2)

mt_plots_PCA <- function(
  D,          # SummarizedExperiment input
  title="",   # title of plot
  ...         # additional arguments directly passed to aes() of ggplot
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # extract data and verify that there are no NA values
  X = t(assay(D))
  if (any(is.na(X))) stop("Data matrix for PCA cannot contain NAs")
  
  # PCA
  pca <- prcomp(x=as.matrix(X)) 
  expvar = (pca$sdev)^2 / sum(pca$sdev^2) 
  # assemble data frame, two PCS and sample info
  df = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], colData(D))
  
  # prep
  pc1name = sprintf('PC1 (%.1f%%)', expvar[1]*100)
  pc2name = sprintf('PC2 (%.1f%%)', expvar[2]*100)
  
  # plot
  p <- ggplot(data=df) +
    geom_point(aes(x=PC1,y=PC2,...)) + 
    xlab(pc1name) + ylab(pc2name) + ggtitle(title)
  
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("PCA, aes: %s", mti_dots_to_str(...)),
      output = p
    )
  
  # return
  D
  

}
