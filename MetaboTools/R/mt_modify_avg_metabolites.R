#' Averages duplicate metabolites
#'
#' Averages values of metabolites (rows) with same values in specified rowData column.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param group_col Name of rowData column (metabolite annotation) by which duplicates can be identified.
#'
#' @return assay: duplicate metabolites (rows) combined
#'
#' @examples
#' \dontrun{... %>% mt_modify_avg_metabolites(group_col = 'name') %>% ...}
#'
#' @author RB
#'
#' @import dplyr
#'
#' @export
mt_modify_avg_metabolites <- function(D, group_col) {

  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(group_col))
    stop("group_col can't be empty")

  total_rows <- nrow(D)
  # TODO: check that rowdata for duplicates are the same, throw mti_logwarning if not

  # mean per metabolite by group_col column in rowData
  X <- dplyr::bind_cols(D %>% assay() %>% data.frame(), D %>% rowData() %>% data.frame() %>%
                          dplyr::select(!!as.name(group_col))) %>%
    dplyr::group_by(!!as.name(group_col)) %>%
    dplyr::summarise_at(dplyr::vars(-group_cols()), mean) %>% # mean per metabolite
    dplyr::ungroup()

  # identify the unique rows
  unique_rows <- D %>% rowData() %>% data.frame() %>% dplyr::group_by(!!as.name(group_col)) %>%
    group_rows() %>%  lapply(., FUN=function(x)x[1]) %>% unlist()

  # subset SE for unique rows
  D <- D[unique_rows, ]

  # replace assay data with means
  assay(D, withDimnames = F) <- X[match(rowData(D)[[group_col]], X[[group_col]]), ] %>% select(-!!as.name(group_col))

  ## add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('%d duplicate metabolites averaged', (total_rows-unique_rows))
    )

  ##return
  D

}
