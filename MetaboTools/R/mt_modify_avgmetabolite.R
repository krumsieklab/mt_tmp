#' mt-modify_avgmetabolite
#'
#' Averages duplicate metabolites
#'
#' @param D \code{SummarizedExperiment} input
#' @param avg_by name of rowData column (metabolite annotation) by which duplicates can be identified
#'
#' @return D with duplicate metabolites combined
#'
#' @examples
#' \dontrun{... %>% mt_modify_avgmetabolite(avg_by = 'name') %>% ...}
#'
#' @author RB
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_modify_avgmetabolite <- function(
  D,       # SummarizedExperiment input
  avg_by   # metabolite annotation column to compare with
) {

  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(avg_by))
    stop("avg_by can't be empty")

  total_rows <- nrow(D)
  # TODO: check that rowdata for duplicates are the same, throw mti_logwarning if not
  
  # mean per metabolite by avg_by column in rowData
  X <- dplyr::bind_cols(D %>% assay() %>% data.frame(), D %>% rowData() %>% data.frame() %>% 
                          select(!!as.name(avg_by))) %>% 
    dplyr::group_by(!!as.name(avg_by)) %>% 
    dplyr::summarise_at(vars(-group_cols()), mean) %>% # mean per metabolite 
    ungroup()
  
  # identify the unique rows
  unique_rows <- D %>% rowData() %>% data.frame() %>% dplyr::group_by(!!as.name(avg_by)) %>% 
    group_rows() %>%  lapply(., FUN=function(x)x[1]) %>% unlist()
 
  # subset SE for unique rows
  D <- D[unique_rows, ]
  
  # replace assay data with means
  assay(D) <- X[match(rowData(D)[[avg_by]], X[[avg_by]]), ] %>% select(-!!as.name(avg_by))
  
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
