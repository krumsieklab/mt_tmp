#' Create new variable by expression
#'
#' Creates a new variable using dplyr's mutate() mechanism. Formula can use any field in the respective colData/rowData
#'
#' @param D \code{SummarizedExperiment} input
#' @param anno_type direction, "samples" or "metabolites"
#' @param col_name name of new variable
#' @param term mutate term to forward to dplyr::mutate
#'
#' @return SummarizedExperiment with altered colData or rowData, containing a new field.
#'
#' @examples
#' \dontrun{# Convert numeric sample annotation field 'num' to factor
#' ...  %>%
#'  mt_modify_mutate(anno_type="samples", col_name="num", term = as.factor(num)) %>% ...}
#'
#' @author JK
#'
#' @importFrom data.table :=
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_modify_mutate <- function(D, anno_type, col_name, term) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(anno_type %in% c("samples","metabolites"))) stop("anno_type must be either 'samples' or 'metabolites'")

  # run mutate argument in correct direction
  x <- dplyr::enquo(term)

  if (anno_type=="samples") {
    cn <- colnames(D) # mutate destroys colnames
    cd <- colData(D) %>% as.data.frame()
    cd %<>% dplyr::mutate(!!col_name := !!x)
    colData(D) <- DataFrame(cd)
    colnames(D) <- cn
  } else if (anno_type=="metabolites") {
    cn <- colnames(D) # mutate destroys colnames
    rd <- rowData(D) %>% as.data.frame()
    rd %<>% dplyr::mutate(!!col_name := !!x)
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
      logtxt = sprintf("added variable to %s: %s := %s", anno_type, col_name, as.character(x))
    )

  ## return
  D

}
