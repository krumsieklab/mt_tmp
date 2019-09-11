#' Output information about dataset to log.
#'
#' @param D 
#'
#' @return Nothing - leaves \code{D} unchanged
#'
#' @examples
#'...  %>% 
#'   mt_logging_tic() %>% ...
#' 
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