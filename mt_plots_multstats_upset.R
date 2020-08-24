#' Upset plot of Multiple Statistical Objects
#'
#' Generate an upset overlap plot of the significant results of multiple statistical objects.
#'
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_list List of stat names to plot; default NA (for all stat objects in D)
#' @param cutoff p.adj significance cutoff
#'
#' @return results$output: upset plot of overlap
#'
#'#' @examples
#' \dontrun{D <- D %>% mt_plots_pheatmap_multstats(cutoff=0.05)}
#'
#' @author JK
#'
#' @export

mt_plots_multstats_upset <- function(
  D, stat_list=NA, cutoff) {

  # result collection part of this code is copied from KC's mt_plots_multstats_heatmap

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # if NA, stat_list is all stat results
  if(is.na(stat_list)){
    stat_list <- MetaboTools:::mti_res_get_stats_entries(D)%>%
      purrr::map_chr(~.x$output$name)
  }

  # verify that there are results
  if(length(stat_list) == 0){
    stop("There are no statistical results.")
  }

  # extract stat results
  stat_res <- stat_list %>% lapply(function(i){MetaboTools:::mti_get_stat_by_name(D, i) %>%
      dplyr::mutate(stat_name = i)
  })

  # generate list of significant metabolites
  signif_list <- stat_res %>% lapply(function(x){
     x %>% filter(p.adj<cutoff) %>% .$var
  })
  names(signif_list) <- stat_list

  # upset plot
  # UpSetR::upset(UpSetR::fromList(signif_list), nsets=length(signif_list), decreasing=c(T,T), order.by=c("degree", "freq"))

  load("signif_list.rds")
  signif_list %>% ComplexHeatmap::make_comb_mat(mode="distinct") %>% ComplexHeatmap::UpSet()
  signif_list %>% ComplexHeatmap::make_comb_mat(mode="intersect") %>% ComplexHeatmap::UpSet()
}

