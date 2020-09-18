#' Add Missingness Column to rowData / colData
#'
#' Adds a missingness column for rowData (metabolites) and / or colData (samples).
#'
#' @param D \code{SummarizedExperiment} input
#' @param add_mets add metabolite missingness? default: T
#' @param add_samples add sample missingness? default: T
#' @param col_name_mets name to give the new rowData column; default: missing
#' @param col_name_samples name to give the new colData column; default: missing
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
                                  add_samples=T,
                                  col_name_mets="missing",
                                  col_name_samples="missing"){

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
    if(col_name_mets %in% colnames(rowData(D))){
      warning(paste0(col_name_mets, " already exists in rowData. It will be overwritten."))
    }
    miss=missingness(X)
    rowData(D)[[col_name_mets]] <- miss
  }

  # add sample missingness
  if(add_samples){
    if(col_name_samples %in% colnames(colData(D))){
      warning(paste0(col_name_samples, " already exists in colData. It will be overwritten."))
    }
    miss=missingness(t(X))
    colData(D)[[col_name_samples]] <- miss
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
