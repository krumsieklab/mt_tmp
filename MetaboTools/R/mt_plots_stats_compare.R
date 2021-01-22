#' Comparative plot between two comparisons.
#'
#' Produces a plot that compares the directed -log10 p-values between two previously executed stats.
#'
#' @param D1 first SE dataset to compare; the one in the pipeline
#' @param stat1 name of statistical comparison in first dataset
#' @param filter1 filter term, defining which metabolites to label from first comparison (can use elements of stats table)
#' @param D2 second SE dataset, if not given, will be the same as the first
#' @param stat2 name of statistical comparison in second dataset
#' @param filter2 filter term, defining which metabolites to label from second comparison (can use elements of stats table)
#' @param filter_op  if AND -> two colors, one for those where both stats match the criterion, and one where they don't
#'                  if OR -> three colors, a third one where only one stat matches the criterion
#' @param plot_title optional param for plot title
#' @param label_col optional argument on which column in the statistical results df to use for labeling points
#' @param point_size size of the points on the ggplot
#' @param return_plot_only return only the plot object. WARNING: setting this to true makes the function non-MT pipeline compatible.
#' @param export_file WHAT DOES THIS DO?
#' @param use_estimate use estimate for comparison, instead of statistic; default: F
#'
#' @return $result: plot, p-value histogram
#'
#' @examples
#' \dontrun{## compare two stats from inside the same pipeline
#' ... %>%
#' mt_plots_stats_compare(stat1='WT',
#'   filter1= p.adj<0.1,
#'   stat2='KO',
#'   filter2= p.adj<0.1,
#'   filter_op = 'OR'
#' ) %>% ...
#'
#' ## compare two stats from different pipelines, as part of the pipeline of the second
#' # 'comp' is a string that contains the name of a comparison (here both SEs have the same comparison on two datasets)
#' .. %>% mt_plots_stats_compare(
#'   stat1 = comp, filter1 = p.adj<0.1,
#'   D2 = firstPipeSE, stat2 = comp, filter2 = p.adj<0.1,
#'   filter_op = "OR") %>% ...
#'
#' ## compare two stats from different pipelines, output as plot object
#' ## not part of the actual MT pipelines, but separate call
#' # 'comp' is a string that contains the name of a comparison (here both SEs have the same comparison on two datasets)
#' gg <- mt_plots_stats_compare(
#'   D1 = D1, stat1 = comp, filter1 = p.adj<0.1,
#'   D2 = D2, stat2 = comp, filter2 = p.adj<0.1,
#'   filter_op = "OR", return_plot_only=T)
#' }
#'
#' @author JK
#'
#' @import ggplot2
#'
#' @export
mt_plots_stats_compare <- function(D1,
                                   stat1,
                                   filter1,
                                   D2 = D1,
                                   stat2,
                                   filter2,
                                   filter_op="AND",
                                   plot_title = "",
                                   label_col = "name",
                                   point_size = 1.5,
                                   return_plot_only=F,
                                   export_file = NULL,
                                   use_estimate = F) {

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D1))
  stopifnot("SummarizedExperiment" %in% class(D2))
  if (missing(filter1)) stop("Must provide 'filter1'")
  if (missing(filter2)) stop("Must provide 'filter2'")

  filter1q <- dplyr::enquo(filter1)
  filter2q <- dplyr::enquo(filter2)

  if (!(filter_op %in% c("AND","OR"))) stop("filter_op must be 'AND' or 'OR'")

  ## obtain the two stats structures
  s1 <- MetaboTools:::mtm_get_stat_by_name(D1, stat1, fullstruct=T)
  s1t <- s1$table
  s2 <- MetaboTools:::mtm_get_stat_by_name(D2, stat2, fullstruct=T)
  s2t <- s2$table

  # if use_estimate==T, check that estimate column exists
  if(use_estimate){
    if("estimate" %in% colnames(s1t) == F){stop("Column estimate was not found in stat table 1")}
    if("estimate" %in% colnames(s2t) == F){stop("Column estimate was not found in stat table 2")}
  }

  ## add directed p-value, filter, and merge
  s1t$dp1 <- ifelse(use_estimate, -log10(s1t$p.value) * sign(s1t$estimate), -log10(s1t$p.value) * sign(s1t$statistic))
  s1t$filtered1 <- s1t$var %in% (s1t %>% dplyr::filter(!!filter1q))$var
  s2t$dp2 <- ifelse(use_estimate, -log10(s2t$p.value) * sign(s2t$estimate), -log10(s2t$p.value) * sign(s2t$statistic))
  s2t$filtered2 <- s2t$var %in% (s2t %>% dplyr::filter(!!filter2q))$var
  st <- merge(s1t, s2t, by='var')
  st <- merge(st, rowData(D1), by.x="var", by.y='row.names', all.x=T) # add names


  # combine filters
  if (filter_op=="AND")
    # AND
    st$filtered = as.numeric(st$filtered1 & st$filtered2)
  else if (filter_op=="OR")
    # OR
    st$filtered = as.numeric(st$filtered1) + as.numeric(st$filtered2)
  else
    stop("bug")

  ## create axis labels
  if ("groups" %in% names(s1) && length(s1$groups)==2) {
    xlabel = sprintf("%s high <--   dir. log10(p-value)   --> %s high", s1$groups[1], s1$groups[2])
  } else {
    xlabel = 'directed log10(p-value)'
  }
  if ("groups" %in% names(s2) && length(s2$groups)==2) {
    ylabel = sprintf("%s high <--   dir. log10(p-value)   --> %s high", s2$groups[1], s2$groups[2])
  } else {
    ylabel = 'directed log10(p-value)'
  }


  ## plot
  st <- as.data.frame(st)
  p <- st %>%
    ggplot2::ggplot(ggplot2::aes(x=dp1,y=dp2,color=as.factor(filtered))) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(color='filtered') +
    ggrepel::geom_text_repel(data=dplyr::filter(st, filtered>0), ggplot2::aes_string(label=label_col), size=3, colour = "black") +
    ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)

  if (plot_title != "") {
    p <- p + ggplot2::ggtitle(plot_title)

  }

  ## export to file?
  if (!is.null(export_file)) {
    # can't handle list columns, drop those
    keep = !sapply(st, is.list)
    openxlsx::write.xlsx(x=st[,keep], file=export_file, asTable=F)
  }

  if (!return_plot_only) {
    ## add status information & plot
    funargs <- MetaboTools:::mti_funargs()
    metadata(D1)$results %<>%
      MetaboTools::: mti_generate_result(
        funargs = funargs,
        logtxt = sprintf("comparison plot between '%s' and '%s'", stat1, stat2),
        output = list(p),
        output2 = st
      )
    ## return
    D1
  } else {
    p
  }

}
