#' Store heading that will be added to report later.
#' 
#' Will be used when calling \code{mt_reporting_generateMD}.
#'
#' @param D  \code{SummarizedExperiment} input
#' @param text Heading text
#'
#' @return $result: stores info about heading
#'
#' @examples
#' ... %>% 
#' # add first and second level of heading
#' mt_reporting_text("Preprocessing") %>%
#' mt_reporting_text("Data were first filtered for missing values, then normalized and log-transformed") %>%
#' ...
#' 
#' @author EB
#' 
#' @export
mt_reporting_text <- function(
  D,
  text
) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  # validate input
  stopifnot("character" %in% class(text)) 
  
  # add status information & heading info
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = glue("reporting text"),
      output = list(text=text)
    )
  
  # return
  D
  
}


