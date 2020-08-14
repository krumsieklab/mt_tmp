#' Validate MD5 checksum of file.
#'
#' Note: Has to be run after loading the data, since it needs to operate on an existing SummarizedExperiment.
#'
#' @param D SummarizedExperiment for pass-through
#' @param file File path... does not necessarily have to be the same as the dataset loaded in the pipeline
#' @param checksum Checksum to test for.
#'
#' @return Does not change SummarizedExperiment, only adds a log message.
#'
#' @examples
#' # first call, to get the checksum (will crash, deliberately)
#' \dontrun{... %>% mt_files_checksum(file="input.xlsx", checksum="") %>% ...}
#'
#' # copy-paste the correct ('actual') checksum from the error message into the call:
#' \dontrun{... %>% mt_files_checksum(file="input.xlsx", checksum="688048bd1eb9c771be0eb49548a6f947") %>% ...}
#'
#' @author JK
#'
#' @import SummarizedExperiment
#' @importFrom magrittr %<>%
#'
#' @export
mt_files_checksum <- function(D, file, checksum) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # throw error if file does not exist
  if(!file.exists(file)) stop(sprintf("File does not exist: %s", file))

  # calculate checksum
  md5 <- tools::md5sum(file)[[1]]
  # crash if wrong
  if (!missing(checksum)) {
    if (checksum != md5) {
      stop(sprintf("Wrong checksum for %s, expected: %s, actual: %s\n", file, checksum, md5))
    }
    logtxt <- sprintf("Correct checksum for %s: %s", file, md5)
  } else {
    # no checksum given by user, just show checksum of file
    logtxt <- sprintf("Checksum for %s: %s", file, md5)
  }

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = logtxt
    )

  # return
  D


}
