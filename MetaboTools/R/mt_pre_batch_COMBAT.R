#' ComBat batch correction
#'
#' Performs batch correction via ComBat method from \code{sva} package.
#'
#' @param D \code{SummarizedExperiment} input
#' @param batches sample annotation (colData) column name that contains batch assignment
#'
#' @return assay: batch-corrected version
#'
#' @examples
#' \dontrun{... %>% mt_pre_batch_ComBat(batches="BATCH") %>% ...
#' }
#'
#' @author JK
#'
#' @export
mt_pre_batch_ComBat = function(D, batches) {

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = assay(D) # combat uses data probe s sample

  # crash if there are any missing values
  if (any(is.na(X))) {
    stop("ComBat batch correction cannot work with missing/NA values.")
  }

  # get variable
  if (!(batches %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", batches))
  b = colData(D)[[batches]]
  # ensure it's a factor or character vector
  if (!(is.character(b) || is.factor(b))) stop(sprintf("'%s' has to be character, factor, or numeric.", batches))
  b = as.factor(b)

  # run ComBat
  # invisible (capture.output(     # seems to be the only way to silence ComBat() output to stdout and stderr
  #   invisible( capture.output(
      Y <- sva::ComBat(
        dat=X,
        batch=b,
        mod=NULL,
        par.prior=TRUE,
        prior.plots=FALSE)
  #     ), type = c("output")
  #   ), type = c("message")
  # )
  # )
  # )

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("ComBat batch correction, variable: '%s'",batches)
    )

  # return
  assay(D) = Y
  D


}
