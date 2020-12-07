#' Remove metadata results from SummarizedExperiment
#'
#' Removes either the entire result list or all plots from the metadata of a SummarizedExperiment.
#' Can be used to obtain a more lightweight object for further processing.
#'
#' @description
#' strip=="all" will only leave a single result in the object, the one from this function.
#'
#' @description
#' strip=="plots" will set all ggplot objects to NULL, still allowing HTML reports to be generated.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param strip What to remove. Has to be "all" or "plots". Default: "all".
#'
#' @return $results: Either all plots removed or all results entries removed.
#'
#' @examples
#' \dontrun{# at the end of a pipeline
#' ... %>% mt_other_strip_results(strip="plots")
#' ...}
#'
#' @author JK
#'
#' @export
mt_other_strip_results <- function(D, strip="all") {
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(strip %in% c("all", "plots"))

  # delete all plots?
  if (strip=="plots") {
    # indices to plot results
    inds <- D %>% metadata() %>% .$results %>% purrr::map("fun") %>% purrr::map(~"plots" %in% .) %>% unlist() %>% which()
    # delete plots
    for (ind in inds) {
      metadata(D)$results[[ind]]['output'] <- list(NULL) # [] and list() needs to be used, otherwise $output is deleted
    }
  }

  # add status information & plot
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Stripped metadata results: %s", strip)
    )

  # remove everything but the last entry?
  if (strip=="all") {
    l <- D %>% metadata() %>% .$results %>% length()
    metadata(D)$results <-  metadata(D)$results[l]
  }

  # return
  D


}




