#' Box Plot of samples
#'
#' Create a box plot of the samples. Can be colored by factor.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param title title of box plot. Default: "Sample boxplot".
#' @param legend show legend? Default: T.
#' @param ylabel y axis label. Default: "Metabolite concentrations".
#' @param plot_logged Show plot logged? Note: plot will still be logged if data have been logged before in pipeline. Default: F.
#' @param ggadd Further elements/functions to add (+) to the ggplot object. Default: NULL.
#' @param ...  Additional arguments directly passed to aes() of ggplot.
#'
#' @return $results$output: plot, sample box plot
#'
#' @examples
#' \dontrun{## sample boxplot, color by colData 'group' variable, with specific title, on log scale,
#' ... %>% mt_plots_sample_box(color=group, title='after quotient normalization', plot_logged=T) %>% ...
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_sample_box <- function(D,
                                title="Sample boxplot",
                                legend=T,
                                ylabel = "Metabolite concentrations",
                                plot_logged=F,
                                ggadd=NULL,
                                ...) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # plot_logged?
  Dplot = D
  if (plot_logged) {
    assay(Dplot) <- log2(assay(Dplot))
    ylabel = sprintf("%s [log2]", ylabel)
  }

  # merge with sample annotations, only keep the ones that were actually used
  cd <- Dplot %>% colData() %>% as.data.frame() %>% tibble::rownames_to_column("merge.primary")
  keep <- c(mti_extract_variables(quos(...)), "merge.primary")
  cd <- cd[,colnames(cd) %in% keep,drop=F]
  df <- cbind(cd, t(assay(Dplot)))
  # generate ggplot
  p <- df %>%  tidyr::gather(metab, value, dplyr::one_of(rownames(Dplot))) %>%
    ggplot(aes(x = merge.primary, y = value, ...)) +
    geom_boxplot() +
    ylab(ylabel) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # todo add rowname

  # remove legend?
  if (!legend) p = p + theme(legend.position="none")

  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("sample boxplot, aes: %s", mti_dots_to_str(...)),
      output = list(p)
    )

  # return
  D


}
