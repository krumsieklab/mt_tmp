#' Extract results from a "path" of entries in metadata
#'
#' Extracts metadata()$results entries in a given namespace of arbitrary depth. See examples below.
#'
#' @param D \code{SummarizedExperiment} input
#' @param path path, list of strings
#'
#' @returns SummarizedExperiment with corrected data
#'
#' @examples
#' \dontrun{# get all result entries for plots
#' D %>% mtm_res_get_path("plots")
#' # get all result entries for stats
#' D %>% mtm_res_get_path("stats")
#' # check if quotient normalization has been performed
#' q <- D%>% mtm_res_get_path(c("pre","norm","quot"))
#' if (length(q)==0) stop("No quotient normalization performed.")
#' }
#'
#' @author JK
#' @export
mtm_res_get_path <- function(D,path) {
  # [sorry for the code, but it works]s
  labels <- lapply(metadata(D)$results,'[[', 'fun')
  m <- rep(T,length(labels))
  # exclude non-matches
  for (i in 1:length(path)) {
    m <- m & sapply(labels, function(label){
      if (length(label)<i) F
      else {
        if (label[i]!=path[i])F
        else T
      }
    })
  }
  # return
  metadata(D)$results[m]
}

#' Extract all "pre" $results entries from pipeline
#'
#' @param D SummarizedExperiment
#'
#' @noMd
mtm_res_get_pre_entries  <- function(D){mtm_res_get_path(D,"pre")}

#' Extract all "plots" $results entries from pipeline
#'
#' @param D SummarizedExperiment
#'
#' @noMd
mtm_res_get_plots_entries <- function(D){mtm_res_get_path(D,"plots")}

#' Extract all stats entries.
#'
#' Returns all stats entries from the $results list.
#'
#' @param D \code{SummarizedExperiment} input
#'
#' @return list if stats results
#'
#' @author JK
#' @export
mtm_res_get_stats_entries <- function(D){mtm_res_get_path(D,"stats")}

#' Return stats output by name
#'
#' Finds a named statistical result from a mt_stats... function and returns the $output$table dataframe.
#'
#' @param D SummarizedExperiment
#' @param name Name of statistical comparison
#' @param fullstruct optional, output entire $output structure, not just the $table inside.
#'
#' @return $output$table dataframe
#' @export
mtm_get_stat_by_name <- function(D, name, fullstruct=F){
  stopifnot("SummarizedExperiment" %in% class(D))

  if(! ("results" %in% names(metadata(D))))
    stop("no results element found in D")

  stats <- mtm_res_get_stats_entries(D)

  if(length(stats) == 0)
    stop("no stats element found in D")

  names  <- stats %>% purrr::map_chr(~.x$output$name)
  output <- which(names == name)

  if(length(output) == 0)
    stop("stat element with name ", name, " does not exist")
  if(length(output)  > 1)
    stop("there are multiple stat elements with name ", name)

  if(!fullstruct) {
    output <- stats[[ output ]]$output$table
    if( ! any(c("tibble", "data.frame") %in% class(output)) )
      stop("output of stat result ", stat, " is not a table")

    if( ! ("var" %in% colnames(output)) )
      stop("output of stat result ", name, " does not have 'var' column")
  } else {
    output <- stats[[ output ]]$output
  }
  output
}

#' Extract all plot objects from pipeline
#'
#' Differs from mtm_res_get_plots_entries in that that function extracts full result structures, and this one returns just the plots
#'
#' @param D SummarizedExperiment
#' @param unlist Unlist all plots into one long list (default: T)
#'
#' @noMd
mtm_res_get_plots <- function(D,unlist=T){
  l=sapply(mtm_res_get_path(D,"plots"),'[[','output',simplify=F)
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

