#' 2D PCA of samples.
#'
#' Can be colored by any variable in colData.
#'
#' @param D \code{SummarizedExperiment} input
#' @param title Title of the plot, default: "PCA"
#' @param scaledata  scale data before plotting? (mean 0, std 1), default: FALSE
#' @param PCa first PC to plot, default is 1 (PC1)
#' @param PCb second PC to plot, default is 2 (PC2)
#' @param show show 'scores' or 'loadings'
#' @param labelby field to label. default: none
#' @param textrepel try to avoid all text overlaps when labeling? default:T
#' @param ellipse confidence interval for ellipse. default: none (no ellipse)
#' @param expvarplot add explained variance plot? default: F
#' @param store.matrices store scores and loadings matrices in result structure? default: F
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... # additional expression directly passed to aes() of ggplot, can refer to colData
#'
#' @return result output: plot(s)
#' @return result output2: scores and loadings matrix
#
#' @examples
#' \dontrun{## PCA on scaledata, color and shape by "Group" variable in colData
#' ... $>$ mt_plots_PCA(scaledata=T, color=Group, shape=Group, title="PCA - scaled data") %>% ...
#' ## PCA scores plot on non-scaled data, with ellipse and extra explained variance plot, and two ggadds (white background and centering of title)
#' mt_plots_PCA(title="PCA scores", show = 'scores', scaledata=F, PCa=1, PCb=2, ellipse=0.95, expvarplot=T, ggadd = theme_bw() + theme(plot.title=element_text(hjust=0.5)))
#' }
#'
#' @author JK, BGP
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_PCA <- function(
  D,
  title="PCA",
  scaledata=F,
  PCa=1,
  PCb=2,
  show='scores',
  labelby='',
  textrepel=T,
  ellipse=NA,
  expvarplot=F,
  store.matrices=F,
  ggadd=NULL,
  ...
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # helper function to combine two aesthetics, e.g. from aes() and aes_string()
  combine_aes <- function(...) {
    v <- c(...)
    class(v) <- "uneval"
    v
  }

  # extract data and verify that there are no NA values
  X = t(assay(D))
  if (any(is.na(X))) stop("Data matrix for PCA cannot contain NAs")
  # check that show is either "scores" or "loadings"
  if (!(show %in% c("scores","loadings"))) stop("Show must be either 'scores' or 'loadings'")

  # scale?
  if (scaledata) X <- scale(X) #By default, the scale R-function: mean-centers and scales to unit variance the X matrix

  # PCA
  pca <- stats::prcomp(x=as.matrix(X), center=F, scale=F)
  expvar = (pca$sdev)^2 / sum(pca$sdev^2)
  # assemble data frame, two PCS and sample info
  if (show=='scores') df = data.frame(x = pca$x[,PCa], y = pca$x[,PCb], colData(D)) # scores and colData
  else  df = data.frame(x = pca$rotation[,PCa], y = pca$rotation[,PCb], rowData(D)) # loadings and rowData
  colnames(df)[1:2] <- c(sprintf("PC%d",PCa),sprintf("PC%d",PCb))

  # prep
  pc1name = sprintf('PC%d (%.1f%%)', PCa, expvar[PCa]*100)
  pc2name = sprintf('PC%d (%.1f%%)', PCb, expvar[PCb]*100)

  # plot
  p <- ggplot(data=df) +
    geom_point(combine_aes(aes_string(x=sprintf("PC%d",PCa),y=sprintf("PC%d",PCb)),aes(...))) +
    xlab(pc1name) + ylab(pc2name) + ggtitle(title)
  # add text?
  if (nchar(labelby)>0) {
    if (textrepel) p <- p + ggrepel::geom_text_repel(aes_string(x=sprintf("PC%d",PCa),y=sprintf("PC%d",PCb),label=labelby,...))
    else p <- p + geom_text(combine_aes(aes_string(x=sprintf("PC%d",PCa),y=sprintf("PC%d",PCb),label=labelby),aes(...)))
  }
  # add ellipse?
  if (!is.na(ellipse)) p <- p + stat_ellipse(aes_string(x=sprintf("PC%d",PCa),y=sprintf("PC%d",PCb)), level=ellipse)
  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  # add explained variance plot?
  plotlist <- list(p)
  if (expvarplot) {
    nlimit <- 10
    expdf <- data.frame(n=1:length(expvar), expvar=expvar*100)
    # cut to first nlimit at most
    expdf <- expdf[1:min(nlimit,nrow(expdf)),]
    newp <- ggplot(data=expdf, aes(x=n, y=expvar)) +
      geom_bar(stat="identity", fill="steelblue")+
      ylab("") + xlab("PC") + ggtitle(sprintf("Explained variance (first %d PCs)", nrow(expdf))) +
      geom_text(aes(label=sprintf('%.1f%%',expdf$expvar)), vjust=-0.3, size=3.5) +
      scale_x_discrete(limits=1:nrow(expdf))
    # add to list of plots for results
    plotlist[[2]] <- newp
  }

  # prep output matrices
  if (store.matrices) {
    scores=pca$x
    loadings=pca$rotation
    rownames(loadings) <- D %>% rowData %>% .$name
    output2 <- list(scores=scores, loadings=loadings)
  } else {
    output2 <- NULL
  }

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("PCA, %s, labelby: %s, aes: %s", show, labelby,  mti_dots_to_str(...)),
      output = plotlist,
      output2 = output2
    )

  # return
  D


}
