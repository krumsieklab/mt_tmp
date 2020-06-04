#' Filter metabolites.
#'
#' Filters metabolites according to an expression. Expression can access entries of rowData.
#'
#' @param D \code{SummarizedExperiment} input
#' @param metab_filter Logical expression. Can use fields from \code{rowData()}.
#'
#' @return assay: filtered data matrix with removed metabolites
#' @return rowData: filtered down
#'
#' @examples
#' \dontrun{... %>% mt_modify_filter_metabolites(SUPER_PATHWAY=="Nucleotide") %>% ...}
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_modify_filter_metabolites <- function(D, metab_filter){

    stopifnot("SummarizedExperiment" %in% class(D))
    if(missing(metab_filter))
        stop("metabolite filter can't be empty")

    ## APPLY FILTER TO ROW DATA
    metab_filter_q <- dplyr::enquo(metab_filter)
    rd <- rowData(D) %>%
        as.data.frame() %>%
        dplyr::mutate(rownames = rownames(D)) %>%
        dplyr::filter(!!metab_filter_q)

    ## SUBSET METABOLITES
    excluded <- rownames(D)[ !(rownames(D) %in% rd$rownames) ]
    D <- D[rd$rownames, ]

    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>%
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Filter metabolites: %s",  as.character(metab_filter_q)),
                      output = excluded
                  )
    ## return
    D
}

#' Filter samples.
#'
#' Filters samples according to an expression. Expression can access entries of colData.
#'
#' @param D \code{SummarizedExperiment} input
#' @param sample_filter Logical expression. Can use fields from \code{colData()}.
#'
#' @return assay: filtered data matrix with removed samples
#' @return colData: filtered down
#'
#' @examples
#' # filter to two specific groups of samples
#' \dontrun{... %>% mt_modify_filter_samples(sample_filter = GROUP %in% c("FL","ctrl")) %>% ...}
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @author JK
#'
#' @export
mt_modify_filter_samples <- function(D, sample_filter){

    stopifnot("SummarizedExperiment" %in% class(D))
    if(missing(sample_filter))
        stop("sample filter can't be empty")

    ## APPLY FILTER TO ROW DATA
    sample_filter_q <- dplyr::enquo(sample_filter)
    cd <- colData(D) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("colnames") %>%
        dplyr::filter(!!sample_filter_q)

    ## SUBSET SAMPLES
    cnames <- colData(D) %>% as.data.frame() %>% tibble::rownames_to_column("colnames") %>% .$colnames
    D <- D[, cnames %in% cd$colnames]

    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>%
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Filter samples: %s",  as.character(sample_filter_q)),
                      output = list(kept=cnames %in% cd$colnames)
                  )
    ## return
    D
}
