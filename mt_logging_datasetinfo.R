#' Output information about dataset to log.
#' 
#' Leaves a log entry containing the number of samples, metabolites and annotation fields at the current stage of the pipeline.
#'
#' @param D \code{SummarizedExperiment} input
#'
#' @return Nothing - leaves \code{D} unchanged
#'
#' @examples
#'...  %>% 
#'   mt_logging_datasetinfo() %>% ...
#' 
#' @author JK
#' 
#' @export
mt_logging_datasetinfo <- function(D) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # toc 
  logtxt <- sprintf("Dataset info: %d samples, %d metabolites; %d sample annotation fields, %d metabolite annotation fields",
                    ncol(D), nrow(D), ncol(colData(D)), ncol(rowData(D)))
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = logtxt
    )
  
  # return
  D
  
}