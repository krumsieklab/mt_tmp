#' Wrapper function to filter missingness
#'
#' Filtering missingness of both metabolites and samples, and generating 3 missingness plots: before filtering, after filter metabolites, and after filter samples
#'
#' @param D \code{SummarizedExperiment} input
#' @param plot_options A list of parameters for \code{mt_plots_qc_missingness}
#' @param filter_options A list of oparameters for \code{mt_pre_filtermiss}, default thresholds for metabolites and samples are 0.2 and 0.1
#'
#' @return D with filtered data
#'  column "Met_Missing" will be added to rowData with missing information on metabolites
#'  column "Sample_Missing" will be added to colData with missing information on samples
#' 
#' @examples
#' \dontrun{... %>% mtw_missingness() %>% ...}
#' \dontrun{... %>% mtw_missingness(filter_options = list(met_max = 0.1, sample_max = 0.05)) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

## Add all options of each function used inside
mtw_missingness <- function (D,
                             plot_options = list(),
                             filter_options = list()) {
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # function to merge default and user input
  map_lists <- function(def, user) {
    # check if there are any illegal entries
    bad_args <- setdiff(names(user), names(def))
    if (length(bad_args) > 0) {
      stop(sprintf(
        "Invalid argument(s): %s",
        paste0(bad_args, collapse = ", ")
      ))
    }
    # now fill up entries in def
    def[names(user)] <- user
    # return
    def
  }
  
  plot_options_def = list(
    met_max = NA,
    samp_max = NA,
    plot_mets = T,
    plot_samples = F,
    sec_axis_mets = F,
    sec_axis_samples = F,
    sample_labels = NA,
    plot_data = T
  )
  
  filter_options_def = list(
    met_max = 0.2,
    sample_max = 0.1,
    met_group = NA,
    log_filtered = F,
    log_sample_col = ""
  )
  
  # missingness plot before filtering
  plot_options <- map_lists(plot_options_def, plot_options)
  plot_options$D <- D
  D <- do.call("mt_plots_qc_missingness", plot_options)
  # filtering metabolites?
  filter_options <- map_lists(filter_options_def, filter_options)
  sample_max <- filter_options$sample_max
  filter_options$sample_max <- NULL
  filter_options$D <- D
  D <- do.call("mt_pre_filtermiss", filter_options)
  # plot missingness after filtering metabolites?
  plot_options$D <- D
  D <- do.call("mt_plots_qc_missingness", plot_options)
  # filtering samples?
  filter_options$sample_max <- sample_max
  filter_options$D <- D
  D <- do.call("mt_pre_filtermiss", filter_options)
  # plot missingness after filtering samples?
  plot_options$D <- D
  D <- do.call("mt_plots_qc_missingness", plot_options)
  
  #annotations
  D %<>%
  mt_anno_missingness(anno_type = "metabolites", out_col = "Met_Missing") %>%
  mt_anno_missingness(anno_type = "samples", out_col = "Sample_Missing")
  
  D
}
