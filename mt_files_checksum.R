#' Validate MD5 checksum of file.
#' 
#' Note: Has to be run after loading the data, since it needs to operate on an existing SummarizedExperiment.
#'
#' @param D SummarizedExperiment for pass-through
#' @param file FIle path... does not necessarily have to be the same as the dataset loaded in the pipeline
#' @param checksum Checksum to test for.
#'
#' @return Does not change SummarizedExperiment, only adds a log message.
#'
#' @examples
#' # first call, to get the checksum (will crash, deliberately)
#' ... %>% mt_files_checksum(file="input.xlsx", checksum="") %>% ...   
#' 
#' # copy-paste the correct ('actual') checksum from the error message into the call:
#' ... %>% mt_files_checksum(file="input.xlsx", checksum="688048bd1eb9c771be0eb49548a6f947") %>% ...   
#' 
#' @author JK
#' 
#' @export
mt_files_checksum <- function(D, file, checksum) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # calculate checksum
  md5 <- tools::md5sum(file)[[1]]
  # crash if wrong
  if (checksum != md5) {
    stop(sprintf("Wrong checksum for %s, expected: %s, actual: %s\n", file, checksum, md5))
  } 
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Correct checksum for %s: %s", file, md5)
    )
  
  # return
  D
  
  
}