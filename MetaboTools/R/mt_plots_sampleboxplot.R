#' Boxplot of samples.
#'
#' Can be colored by factor.
#'
#' @param D \code{SummarizedExperiment} input
#' @param plottitle title of boxplot, default "Sample boxplot"
#' @param legend show legend? default: T
#' @param ylabel y axis label, default: "Metabolite concentrations"
#' @param logged show plot logged? default: F. note: plot will still be logged if data have been logged before in pipeline.
#' @param ggadd further elements/functions to add (+) to the ggplot object
#' @param ...  additional arguments directly passed to aes() of ggplot
#'
#' @return $result: plot, sample boxplot
#'
#' @examples
#' \dontrun{## sample boxplot, color by colData 'group' variable, with specific title, on log scale,
#' ... %>% mt_plots_sampleboxplot(color=group, plottitle='after quotient normalization', logged=T) %>% ...
#' }
#'
#' @author JK
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_sampleboxplot <- function(
  D,         # SummarizedExperiment input
  plottitle="Sample boxplot",
  legend=T,  # keep legend?  [could be removed]
  ylabel = "Metabolite concentrations",  # y axis label
  logged=F,  # plot logged
  ggadd=NULL,   # further elements/functions to add (+) to the ggplot object
  ...        # additional arguments directly passed to aes() of ggplot
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # logged?
  Dplot = D
  if (logged) {
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
    ggtitle(plottitle) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # todo add rowname

  # remove legend?
  if (!legend) p = p + theme(legend.position="none")

  # add custom elements?
  if (!is.null(ggadd)) p <- p+ggadd

  # fix ggplot environment
  if (D %>% mti_get_setting("ggplot_fix")) p <- mti_fix_ggplot_env(p)

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
