#' Save SummarizedExperiment to file
#'
#' Simple helper function to save SE to a file. Object name will be \code{D}.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output filename to write to.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @examples
#' \dontrun{#
#' ... %>% mt_write_SE_rds(file="out.RDS") %>% ...
#' ...}
#'
#' @author JK
#'
#' @export
mt_write_SE_rds <- function(D, file) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # write
  save(D, file=file)

  # add status information & plot
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("SummarizedExperiment saved to: %s", file)
    )

  # return
  D

}
