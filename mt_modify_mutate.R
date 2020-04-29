#' Title
#'
#' @param D \code{SummarizedExperiment} input
#' @param dir direction, "samples" or "metabolites"
#' @param varname name of new variable
#' @param term mutate term to forward to dplyr::mutate
#'
#' @return
#' @export
#'
#' @examples
mt_modify_mutate <- function(D, dir, varname, term) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(dir %in% c("samples","metabolites"))) stop("dir must be either 'samples' or 'metabolites'")
  
  # run mutate argument in correct direction
  x <- enquo(term)
  
  if (dir=="samples") {
    cn <- colnames(D) # mutate destroys colnames
    cd <- colData(D) %>% as.data.frame()
    cd %<>% dplyr::mutate(!!varname := !!x)
    colData(D) <- DataFrame(cd)
    colnames(D) <- cn
  } else if (dir=="metabolites") {
    cn <- colnames(D) # mutate destroys colnames
    rd <- rowData(D) %>% as.data.frame()
    rd %<>% dplyr::mutate(!!varname := !!x)
    rowData(D) <- DataFrame(rd)
    colnames(D) <- cn
  } else {
    stop('bug')
  }
  
  
  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("added variable to %s: %s := %s", dir, varname, as.character(x))
    )
  
  ## return
  D
  
}