#' Store heading that will be added to report later.
#'
#' Will be used when calling \code{mt_reporting_generateMD}.
#'
#' @param D  \code{SummarizedExperiment} input (missing if first step in pipeline)
#' @param strtitle Heading text
#' @param lvl Heading level, default: 1 (can be used for nested outline structures)
#'
#' @return $result: stores info about heading
#'
#' @examples
#' \dontrun{... %>%
#' # add first and second level of heading
#' mt_reporting_heading("Preprocessing") %>%
#' mt_reporting_heading("Part 1", lvl=2) %>%
#' ...}
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#' @import ggplot2
#' @export

mt_reporting_heading <- function(
  D,
  strtitle,
  lvl=1
) {

  # if first step in pipeline, create SE
  if(missing(D)){
    # create an empty SummarizedExperiment object
    D <- SummarizedExperiment()
  }else{
    # validate argument
    stopifnot("SummarizedExperiment" %in% class(D))
  }

  # add status information & heading info
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue("reporting heading, level {lvl}: {strtitle}"),
      output = list(lvl=lvl,title=strtitle)
    )

  # return
  D

}


