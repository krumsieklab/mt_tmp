#' 2D UMAP of samples.
#'
#' Can be colored by any variable in colData.
#'
#' @param D \code{SummarizedExperiment} input
#' @param title Title of the plot, default: "UMAP"
#' @param scaledata  scale data before plotting? (mean 0, std 1), default: FALSE
#' @param labelby field to label. default: none
#' @param textrepel try to avoid all text overlaps when labeling? default:T
#' @param store.matrices store scores in result structure? default: F
#' @param n_neighbors ADD DESCRIPTION
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ... # additional expression directly passed to aes() of ggplot, can refer to colData
#'
#' @return result output: plot(s)
#' @return result output2: scores and loadings matrix
#
#' @examples
#' \dontrun{## UMAP on scaledata, color and shape by "Group" variable in colData
#' ... $>$ mt_plots_UMAP(scaledata=T, color=Group, shape=Group, title="UMAP - scaled data") %>% ...
#' }
#'
#' @author JK
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_UMAP <- function(
  D,
  title="UMAP",
  scaledata=F,
  labelby='',
  textrepel=T,
  store.matrices=F,
  ggadd=NULL,
  n_neighbors=15,
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
  if (any(is.na(X))) stop("Data matrix for UMAP cannot contain NAs")

  # scale?
  if (scaledata) X <- scale(X) #By default, the scale R-function: mean-centers and scales to unit variance the X matrix

  # UMAP
  um <- umap::umap(d=as.matrix(X), n_neighbors=n_neighbors)

  # assemble data frame, two components and sample info
  df = data.frame(x = um$layout[,1], y = um$layout[,2], colData(D)) # scores and colData
  colnames(df)[1:2] <- c("comp1","comp2")

  # plot
  p <- ggplot(data=df) +
    geom_point(combine_aes(aes_string(x="comp1",y="comp2"),aes(...))) +
    xlab("comp 1") + ylab("comp 2") + ggtitle(title)
  # add text?
  if (nchar(labelby)>0) {
    if (textrepel) p <- p + ggrepel::geom_text_repel(aes_string(x="comp1",y="comp2",label=labelby,...))
    else p <- p + geom_text(combine_aes(aes_string(x="comp1",y="comp2",label=labelby),aes(...)))
  }

  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  # fix ggplot environment
  if (D %>% mti_get_setting("ggplot_fix")) p <- mti_fix_ggplot_env(p)

  # prep output matrices
  if (store.matrices) {
    output2 <- list(layout=umap$layout)
  } else {
    output2 <- NULL
  }

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("UMAP, labelby: %s, aes: %s", labelby,  mti_dots_to_str(...)),
      output = list(p),
      output2 = output2
    )

  # return
  D


}
