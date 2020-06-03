#' Extract all "pre" $results entries from pipeline
#'
#' @param D SummarizedExperiment
#'
#' @noMd
mtm_res_get_pre_entries  <- function(D){mti_res_get_path(D,"pre")}

#' Extract all "plots" $results entries from pipeline
#'
#' @param D SummarizedExperiment
#'
#' @noMd
mtm_res_get_plots_entries <- function(D){mti_res_get_path(D,"plots")}

#' Extract all "stats" $results entries from pipeline
#'
#' @param D SummarizedExperiment
#'
#' @noMd
mtm_res_get_stats_entries <- function(D){mti_res_get_path(D,"stats")}

#' Extract all plot objects from pipeline
#'
#' Differs from mtm_res_get_plots_entries in that that function extracts full result structures, and this one returns just the plots
#'
#' @param D SummarizedExperiment
#' @param unlist Unlist all plots into one long list (default: T)
#'
#' @noMd
mtm_res_get_plots <- function(D,unlist=T){
  l=sapply(mti_res_get_path(D,"plots"),'[[','output',simplify=F)
  if(unlist) l <- unlist(l,recursive=F)
  l
}

#' Plot all plots from a pipeline
#'
#' Opens a device (default: PDF), plots all plots, closes device.
#' Works either on list of plots, or on SE
#'
#' @param input List of plots or SummarizedExperiment
#' @param dev Device to plot to (default: PDF)
#' @param ... Further paramaters to be passed to dev() function.
#'
#'
#' @noMd
mtm_plot_all_tofile <- function(input, dev=pdf, ...) {
  if ("SummarizedExperiment" %in% class(input)) {
    plots <- mti_multi_get_unique_plots(input)
  } else {
    plots <- input
  }
  dev(...)
  sapply(plots,plot)
  dev.off()
}

