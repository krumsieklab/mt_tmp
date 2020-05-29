#' mt-modify_averagesample
#'
#' Averages duplicate samples.
#'
#' @param D \code{SummarizedExperiment} input
#' @param group_by name of colData column (sample annotation) by which duplicates can be identified
#'
#' @return D with duplicate samples combined
#'
#' @examples
#' \dontrun{... %>% mt_modify_averagesample(group_by = "RID") %>% ...}
#'
#' @author Annalise Schweickart
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_modify_averagesample <- function(
  D,       # SummarizedExperiment input
  group_by   # sample annotation column to compare with
) {

  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(group_by))
    stop("group_by can't be empty")

  # TODO: check that coldata for duplicates are the same, throw mti_logwarning if not
  X <- t(assay(D))
  ave_col <- colData(D)[[group_by]]
  u.col <- sort(unique(ave_col))
  counts <- as.numeric(table(ave_col))
  dup <- unique(ave_col[ave_col%in%u.col[counts > 1]])
  to_remove <- c()
  for (i in 1:length(dup)) {
    index <- which(ave_col == dup[i])
    dat.slice <- X[index, ]
    X[index, ] <- sapply(colMeans(dat.slice, na.rm = TRUE), rep, length(index))
    to_remove <- c(to_remove,index[2:length(index)])
  }
  assay(D) <- t(X)
  D<- D[,-to_remove]

  ## add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf(' %d duplicate samples combined', length(to_remove))
    )

  ##return

  D
}
