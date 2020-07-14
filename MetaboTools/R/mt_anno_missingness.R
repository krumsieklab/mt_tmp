#' Missingness Annotation for Metabolites / Samples
#'
#' Adds a rowData or colData column representing the missingness of metabolites or samples, respectively.
#'
#' @param D \code{SummarizedExperiment} input
#' @param anno_type "samples" (colData) or "metabolites" (rowNames) to add missingness column to
#' @param out_col name of new missingness column (if missing, use "missingness")
#'
#' @return rowData or colData: new annotation column added
#'
#' @examples
#' # example of how to run function
#' \dontrun{
#' #load data
#' mt_files_load_metabolon(file=file, sheet="data") %>%
#' # add metabolite missingness column to rowData
#' mt_anno_missingness(anno_type = "metabolites", out_col = "Metab_Missing")
#' # add sample missingness column to colData
#' mt_anno_missingness(anno_type = "samples", out_col = "Samp_Missing")
#' }
#'
#' @author KC
#'
#' @export
#'

mt_anno_missingness <- function(D,
                                anno_type,
                                out_col){

  # helper function
  # form mt_plots_qc_missingness - should be moved to mt_internal_helpers
  missingness <- function(X)apply(is.na(X),2,sum)/dim(X)[1]

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if (!(anno_type %in% c("samples","metabolites"))) stop("anno_type must be either 'samples' or 'metabolites'")
  if(length(anno_type) > 1) stop("anno_type must be either 'samples' OR 'metabolites' - NOT both")

  # check that assay data hasn't been filtered already?

  X <- t(assay(D))

  if(anno_type == "metabolites"){
    # get missingness of metabolites
    miss_col <- missingness(X)

    # add missingness column to rowData
    if(missing(out_col)){
      out_col <- "missingness"
    }
    rowData(D)[[out_col]] <- miss_col

    logtxt = sprintf("Added missingness annotation for %i metabolites", length(miss_col))

  }else if(anno_type == "samples"){
    # get missingness of samples
    miss_col <- missingness(t(X))

    if(missing(out_col)){
      out_col <- "missingness"
    }
    colData(D)[[out_col]] <- miss_col

    logtxt = sprintf("Added missingness annotation for %i samples", length(miss_col))

  }else{
    stop("Bug")
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
