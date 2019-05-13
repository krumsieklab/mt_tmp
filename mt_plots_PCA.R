require(ggplot2)

#' 2D PCA of samples.
#' 
#' Can be colored by any variable in colData.
#'
#' @param D \code{SummarizedExperiment} input
#' @param title Title of the plot, default: "PCA"
#' @param scaledata  scale data before plotting? (mean 0, std 1), default: FALSE
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... # additional expression directly passed to aes() of ggplot, can refer to colData
#'
#' @return $result: plot, PCA
#
#' @examples
#' ## PCA on scaledata, color and shape by "Group" variable in colData
#' ... $>$ mt_plots_PCA(scaledata=T, color=Group, shape=Group, title="PCA - scaled data") %>% ...
#' 
#' @author JK
#' 
mt_plots_PCA <- function(
  D,           
  title="PCA",    
  scaledata=F,
  ggadd=NULL, 
  ...           
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # extract data and verify that there are no NA values
  X = t(assay(D))
  if (any(is.na(X))) stop("Data matrix for PCA cannot contain NAs")
  
  # scale?
  if (scaledata) X <- scale(X)
  
  # PCA
  pca <- prcomp(x=as.matrix(X), center=F, scale=F) 
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
  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("PCA, aes: %s", mti_dots_to_str(...)),
      output = list(p)
    )
  
  # return
  D
  

}
