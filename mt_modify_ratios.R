#' Generate metabolite ratios.
#' 
#' Transforms the dataset into a new dataset where each 'metabolite' represents a ratio of two metabolites. Warning: For a dataset
#' with originally p metabolites, this will result in p*(p-1) new variables. (e.g. 500 metabolits becomes 249500 ratios).
#' 
#' mt_post_pgain provides a special operation on a ratio data matrix for better interpretation of the resulting p-values.
#'
#' @import SummarizedExperiment
#' @importFrom dplyr %>% mutate gather left_join filter arrange
#'
#' @param D SummarizedExperiment object
#' 
#' @example 
#' # Transform dataset to tratios
#' ... %>%  mt_modify_ratios() %>% ... # proceed with statistical analysis
#'
#' @return SummarizedExperiment containing pairwise ratios from all variables of input
#' @export mt_modify_ratios
#' 
#' @author Jonas Zierer, JK
#' 
mt_modify_ratios <- function(D){

    stopifnot("SummarizedExperiment" %in% class(D))

    as <- assay(D)
    ## as <- matrix(1:(4*2), nrow = 4, ncol = 2, dimnames = list(letters[23:26], letters[1:2]))

    ## FOLDCHANGE FUNCTION (CONSIDER PREVIOUS LOG)
    op <- "/"
    if (length(mti_res_get_path(D, c("pre","trans","log"))) > 0){
        mti_logstatus("data already logscale, using '-'")
        op <- "-"
    }
    
    ## CREATE RATIOS
    as_ratio <- map(1:(nrow(as)-1), ~ sweep(as[(.x+1):nrow(as), , drop = F], 2, as[.x,], "-")) %>%
        setNames(rownames(as)[1:(nrow(as)-1)]) 

    ## CREATE NEW ROWDATA
    rd <- D %>%
        rowData() %>%
        as.data.frame() %>%
        mutate(rownames = rownames(D))
    rd_new <- as_ratio %>%
        imap_dfr(~tibble(m1 = rownames(.x), m2 = .y)) %>%
        mutate(rownames = str_c(m1, m2, sep = "_")) %>%
        left_join(rd %>% transmute(m1 = rownames, name1 = name), by = "m1") %>%
        left_join(rd %>% transmute(m2 = rownames, name2 = name), by = "m2") %>%
        mutate(name = str_c(name1, " / ", name2))
    rd <- rd %>%
        mutate(m1 = rownames,
               m2 = NA,
               name1 = name,
               name2 = NA)
    rd <- bind_rows(rd, rd_new) %>%
        dplyr::select(rownames, m1, m2, name1, name2, everything()) %>%
        column_to_rownames("rownames")
    
    ## COMBINE RATIOS TO SINGLE MATRIX
    as_ratio <- as_ratio %>%
        imap(~{rownames(.x) <- str_c(rownames(.x), .y, sep = "_"); .x}) %>%
        invoke(rbind, .)
    as_ratio <- rbind(as, as_ratio)
    
    ## CHECK NAMES
    if(!identical(rownames(as_ratio), rownames(rd)))
        stop("something went wrong. check data!")

    ## ADD RAW DATA
    rd_final <- bind_rows(rd, rd_new)

    ## CREATE NEW OBJECT
    D <- SummarizedExperiment(assay    = as_ratio,
                              rowData  = rd,
                              colData  = colData(D),
                              metadata = metadata(D))

    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D)$results %<>% 
                  mti_generate_result(
                      funargs = funargs,
                      logtxt = sprintf("Create Metabolite Ratios"),
                      output = NULL
                  )
    D

}






            
