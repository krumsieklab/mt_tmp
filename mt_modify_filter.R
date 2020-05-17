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
#' ... %>% mt_modify_filter_metabolites(SUPER_PATHWAY=="Nucleotide") %>% ...
#' 
#' @author JK
#' 
mt_modify_filter_metabolites <- function(D, metab_filter){

    stopifnot("SummarizedExperiment" %in% class(D))
    if(missing(metab_filter))
        stop("metabolite filter can't be empty")
    
    ## APPLY FILTER TO ROW DATA
    metab_filter_q <- enquo(metab_filter)
    rd <- rowData(D) %>%
        as.data.frame() %>%
        mutate(rownames = rownames(D)) %>%
        filter(!!metab_filter_q)

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
#' Filters metabolites according to an expression. Expression can access entries of colData.
#'
#' @param D \code{SummarizedExperiment} input
#' @param sample_filter Logical expression. Can use fields from \code{colData()}.
#'
#' @return assay: filtered data matrix with removed samples
#' @return rowData: filtered down
#' 
#' @examples
#' # filter to two specific groups of samples
#' ... %>% mt_modify_filter_samples(sample_filter = GROUP %in% c("FL","ctrl")) %>% ...
#' 
#' @author JK
#' 
mt_modify_filter_samples <- function(D, sample_filter){

    stopifnot("SummarizedExperiment" %in% class(D))
    if(missing(sample_filter))
        stop("sample filter can't be empty")
    
    ## APPLY FILTER TO ROW DATA
    sample_filter_q <- enquo(sample_filter)
    cd <- colData(D) %>%
        as.data.frame() %>%
        rownames_to_column("colnames") %>%
        filter(!!sample_filter_q)

    ## SUBSET SAMPLES
    cnames <- colData(D) %>% as.data.frame() %>% rownames_to_column("colnames") %>% .$colnames
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
