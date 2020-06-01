#' Reorder a sample annotation.
#'
#' Does not change anything in the actual dataset, but just the factor order of a column in colData. Can be used to change the plotting order in boxplots, volcano plots, etc.
#'
#' @param D \code{SummarizedExperiment} input
#' @param col_name name of the column in colData to access
#' @param new_order list of strings containing the new order
#'
#' @returns colData: changes the factor order of one column
#'
#' @examples
#' \dontrun{%>% mt_modify_sampleanno_reorder("Group", c('WT_Norm_F','WT_Hyp_F','KO_Norm_F','KO_Hyp_F','WT:KO2_F','WT:WT2_F')) %>%}
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_modify_sampleanno_reorder <- function(
  D,          # SummarizedExperiment input
  col_name,      # name of column to access
  new_order   # new order, either numeric vector with new order, or has to contain all column names
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # get variable
  if (!(col_name %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", col_name))
  p = colData(D)[[col_name]]

  # ensure factor
  p <- as.factor(p)

  # different actions for numeric or factor orders
  if (is.numeric(new_order)) {
    if (length(new_order) != length(levels(p)) || !(all.equal(sort(new_order), 1:length(levels(p)))==T)) stop(sprintf("For numeric new_order, the vector has to contain all numbers from 1 to number of levels. Expected: 1:%d",length(levels(p))))
    p <- gdata::reorder.factor(p, new.order=new_order)
  } else {
    if (length(new_order) != length(levels(p)) || !(all.equal(sort(new_order), sort(levels(p)))==T)) stop(sprintf("For string new_order, the vector has to contain all levels exactly once. Levels: %s", paste0(levels(p),collapse=',')))
    p <- gdata::reorder.factor(p, new.order=new_order)
  }

  # write back
  colData(D)[[col_name]] <- p

  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Reordered column '%s' as '%s'", col_name, paste0(new_order,collapse=','))
    )
  ## return
  D

}






