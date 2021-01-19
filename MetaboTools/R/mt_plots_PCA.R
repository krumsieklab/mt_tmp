#' 2D PCA of samples.
#'
#' Can be colored by any variable in colData.
#'
#' @param D \code{SummarizedExperiment} input
#' @param title Title of the plot, default: "PCA"
#' @param scale_data  scale data before plotting? (mean 0, std 1), default: FALSE
#' @param PC1 first PC to plot, default is 1 (PC1)
#' @param PC2 second PC to plot, default is 2 (PC2)
#' @param data_type Data of type 'scores' or 'loadings'
#' @param label_col field to label. default: none
#' @param text_repel try to avoid all text overlaps when labeling? default:T
#' @param ellipse confidence interval for ellipse. default: none (no ellipse)
#' @param exp_var_plot add explained variance plot? default: F
#' @param store_matrices store scores and loadings matrices in result structure? default: F
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... # additional expression directly passed to aes() of ggplot, can refer to colData
#'
#' @return result output: plot(s)
#' @return result output2: scores and loadings matrix
#
#' @examples
#' \dontrun{## PCA on scale_data, color and shape by "Group" variable in colData
#' ... $>$ mt_plots_pca(scale_data=T, color=Group, shape=Group, title="PCA - scaled data") %>% ...
#' ## PCA scores plot on non-scaled data, with ellipse and extra explained variance plot, and two ggadds (white background and centering of title)
#' mt_plots_pca(title="PCA scores", data_type = 'scores', scale_data=F, PC1=1, PC2=2, ellipse=0.95, exp_var_plot=T, ggadd = theme_bw() + theme(plot.title=element_text(hjust=0.5)))
#' }
#'
#' @author JK, BGP
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_pca <- function(
  D,
  title="PCA",
  scale_data=F,
  PC1=1,
  PC2=2,
  data_type='scores',
  label_col='',
  text_repel=T,
  ellipse=NA,
  exp_var_plot=F,
  store_matrices=F,
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
  # check that data_type is either "scores" or "loadings"
  if (!(data_type %in% c("scores","loadings"))) stop("Show must be either 'scores' or 'loadings'")

  # scale?
  if (scale_data) X <- scale(X) #By default, the scale R-function: mean-centers and scales to unit variance the X matrix

  # PCA
  pca <- stats::prcomp(x=as.matrix(X), center=F, scale=F)
  expvar = (pca$sdev)^2 / sum(pca$sdev^2)
  # assemble data frame, two PCS and sample info
  if (data_type=='scores') df = data.frame(x = pca$x[,PC1], y = pca$x[,PC2], colData(D)) # scores and colData
  else  df = data.frame(x = pca$rotation[,PC1], y = pca$rotation[,PC2], rowData(D)) # loadings and rowData
  colnames(df)[1:2] <- c(sprintf("PC%d",PC1),sprintf("PC%d",PC2))

  # prep
  pc1name = sprintf('PC%d (%.1f%%)', PC1, expvar[PC1]*100)
  pc2name = sprintf('PC%d (%.1f%%)', PC2, expvar[PC2]*100)

  # plot
  p <- ggplot(data=df, combine_aes(aes_string(x=sprintf("PC%d",PC1),y=sprintf("PC%d",PC2)),aes(...))) +
    geom_point() +
    xlab(pc1name) + ylab(pc2name) + ggtitle(title)
  # add text?
  if (nchar(label_col)>0) {
    if (text_repel) p <- p + ggrepel::geom_text_repel(aes_string(x=sprintf("PC%d",PC1),y=sprintf("PC%d",PC2),label=label_col,...))
    else p <- p + geom_text(combine_aes(aes_string(x=sprintf("PC%d",PC1),y=sprintf("PC%d",PC2),label=label_col),aes(...)))
  }
  # add ellipse?
  if (!is.na(ellipse)) p <- p + stat_ellipse(level=ellipse)
  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  # fix ggplot environment
  if (D %>% mti_get_setting("ggplot_fix")) p <- mti_fix_ggplot_env(p)


  # add explained variance plot?
  plotlist <- list(p)
  if (exp_var_plot) {
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
    if (D %>% mti_get_setting("ggplot_fix")) plotlist[[2]] <- mti_fix_ggplot_env(newp)
  }

  # prep output matrices
  if (store_matrices) {
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
      logtxt = sprintf("PCA, %s, label_col: %s, aes: %s", data_type, label_col,  mti_dots_to_str(...)),
      output = plotlist,
      output2 = output2
    )

  # return
  D


}
