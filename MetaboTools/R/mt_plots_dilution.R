#' Correlate variable with dilution factors from quotient normalization.
#'
#' \itemize{
#' \item for factors: will produce boxplot
#' \item for quantitative variable: will produce scatter plot
#' \item requires exactly one previous quotient normalization call in the pipeline
#' }
#'
#' @param D \code{SummarizedExperiment} input
#' @param comp name of colData colum  (sample annotation) to correlation dilution factors with
#' @param boxpl  produce boxplot (TRUE, default), or beeswarm plot (FALSE) .... only relevant for factor comparison
#' @param ggadd further elements/functions to add (+) to the ggplot object
#'
#' @return $result: plot, comparison with dilution factors
#'
#' @examples
#' \dontrun{%>% mt_plots_dilution(comp="group") %>% # compare with 'group' sample annotation
#' }
#'
#' @author JK
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_plots_dilution <- function(
  D,       # SummarizedExperiment input
  comp,    # sample annotation column to compare with
  boxpl=T,  #
  ggadd=NULL
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(length(comp)==1)

  # validate that there has been exactly one quotient normalization step
  q <- D%>% mti_res_get_path(c("pre","norm","quot"))
  if (length(q)>1) stop("There has been more than one quotient normalization call.")
  if (length(q)==0) stop("No quotient normalization performed.")
  # get dilution factors
  vd = q[[1]]$output$dilution

  # get variable to compare to
  if (!(comp %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", comp))
  vc = colData(D)[[comp]]

  # generate data frame for plotting
  dfplot <- data.frame(
    dilution.factor = vd,
    y = vc
  )
  colnames(dfplot)[2] <- comp

  # either produce boxplot or scatter plot
  if (is.character(vc) || is.factor(vc)) {
    # ensure factor
    dfplot[[comp]] = as.factor(dfplot[[comp]])
    # boxplot
    p <- dfplot %>%
      ggplot(aes_string(x=comp,y="dilution.factor",color=comp)) +
      ggtitle("quotient normalization dilution factors")
    if (boxpl)
      p <- p+geom_boxplot()
    if (!boxpl)
      p <- p+ggbeeswarm::geom_quasirandom()

  } else {
    if (!is.numeric(vc)) stop(sprintf("'%s' has to be character, factor, or numeric.", comp))
    # scatter
    p <- dfplot %>% ggplot() +
      geom_point(aes_string(x=comp,y="dilution.factor")) +
      ggtitle("quotient normalization dilution factors")
  }

  if (!is.null(ggadd)) p <- p+ggadd

  # fix ggplot environment
  if (D %>% mti_get_setting("ggplot_fix")) p <- mti_fix_ggplot_env(p)

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("dilution factor plot, '%s'",comp),
      output = list(p)
    )

  # return
  D

}

