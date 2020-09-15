#' Add Missingness Column to rowData / colData
#'
#' Adds a missingness column for rowData (metabolites) and / or colData (samples).
#'
#' @param D \code{SummarizedExperiment} input
#' @param add_mets add metabolite missingness? default: T
#' @param add_samples add sample missingness? default: T
#'
#' @return if add_mets == T: new rowData column
#' @return if add_samples == T: new colData column
#'
#' @examples
#' # example of how to run function
#' \dontrun{}
#'
#' @author KC
#'
#' @export

mt_modify_missingness <- function(D,
                                  add_mets=T,
                                  add_samples=T){

  # helper function
  missingness <- function(X)apply(is.na(X),2,sum)/dim(X)[1]

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  if((add_mets | add_samples) == F){
    stop("At least one of the following boolean argumnets must be TRUE:\nadd_mets\nadd_samples")
  }

  # get data
  X <- t(assay(D))

  # add metabolite missingness
  if(add_mets){
    miss=missingness(X)
    rowData(D)[["missingess"]] <- miss
  }

  # add sample missingness
  if(add_samples){
    miss=missingness(t(X))
    colData(D)[["missingness"]] <- miss
  }

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("missingness columns added, missing values: %d out of %d (%.2f%%)",
                       sum(is.na(X)), nrow(X)*ncol(X), sum(is.na(X))/(nrow(X)*ncol(X))*100)
    )

  # return
  D

}
