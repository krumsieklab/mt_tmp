#' Make sample x pathway abundance table
#'
#' Create a table of sample x pathway abundance
#'
#' @param D \code{SummarizedExperiment} input
#' @param pw_col name of column containing pathways IDs
#' @param size_normalize logical whether to divide summed pathway abundances by
#' size of pathway
#' @param sample_normalize logical whether to divide each sample by its L2
#' norm before processing
#' 
#' @return $result: TODO
#' 
#' @example TODO
#' 
#' @author Massoud Maher
#' 
mt_post_aggregate_pathway <- function(D, pw_col="kegg_db", size_normalize=T,
                                      sample_normalize=T) {

  stopifnot("SummarizedExperiment" %in% class(D))
  meta_D <- metadata(D)
  if(!"pathways" %in% names(meta_D)) {
    stop("'pathways' does not exist in current SummarizedExperiment input")
  }
  if (!pw_col %in% names(meta_D$pathways)) {
    stop(sprintf("'%s' not found in metabolite annotations.", pw_col))
  }
  
  metab.list <- rowData(D) %>% 
    as_tibble(rownames=NA) %>%
    dplyr::select(name, !!rlang::sym(pw_col)) %>%
    unnest(!!rlang::sym(pw_col))
 
  pathway.df <- metab.list %>% 
    group_by(!!rlang::sym(pw_col)) %>%
    dplyr::summarise(pathway.size = n())
  
  n_pathway <- nrow(unique(metab.list[,pw_col]))
  n_sample <- ncol(assay(D))
  
  sample.path <- matrix(nrow=n_sample, ncol=n_pathway)
  samples <- colnames(assay(D))
  pathways <- unlist(unique(metab.list[,pw_col]))
  
  assay.df <- assay(D) %>% as_tibble(rownames="metabolite")
  if(sample_normalize) {
    assay.df <- assay.df %>% mutate_if(is.numeric, function(s) {
      #return(s / norm(s, "2"))
      return(s / sum(s))
    })
  }
  print("assay.df")
  print(assay.df)
  
  for(pi in seq(n_pathway)) {
    path <- pathways[pi]
    path.metabs <- metab.list %>% 
      filter(!!rlang::sym(pw_col) == path) %>%
      dplyr::select(name) %>%
      unlist()
    sample.path.row <- assay.df %>% 
      filter(metabolite %in% path.metabs) %>%
      dplyr::select(-metabolite) %>%
      as.matrix() %>%
      colSums()
    if(size_normalize) {
      size <- pathway.df %>% 
        filter(!!rlang::sym(pw_col) == path) %>% 
        dplyr::select(pathway.size) %>% unlist
      sample.path.row <- sample.path.row / size
    }
    sample.path[,pi] <- sample.path.row
  }
  
  sample.path <- as_tibble(sample.path)
  colnames(sample.path) <- pathways
  sample.path$sample <- samples
  sample.path <- sample.path %>% dplyr::select(sample, everything())
   
   
  return(sample.path)
}

##
#res <- readRDS("/Users/massoudmaher/Documents/Code/fructose/results/mt_no_diff/result_list.Rds")
foo <- mt_post_aggregate_pathway(res$feces)
write_csv(foo, "~/Desktop/norm_feces.csv")
##