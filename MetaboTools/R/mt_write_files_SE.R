#' Save SummarizedExperiment to file
#'
#' Simple helper function to save SE to a file. Object name will be \code{D}
#'
#' @param D SummarizedExperiment input
#' @param file File to write to.
#'
#' @examples
#' \dontrun{#
#' ... %>% mt_write_files_SE(file="out.RDS") %>% ...
#' ...}
#'
#' @author JK
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_write_files_SE <- function(D, file) {

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
