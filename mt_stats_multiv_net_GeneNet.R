#' Computes partial correlation matrix using the GeneNet estimator.
#' 
#' Implementation according to Schäfer and Strimmer, 2006\cr
#' https://www.ncbi.nlm.nih.gov/pubmed/16646851
#'
#' @param D \code{SummarizedExperiment} input
#' @param name name of the correlation matrix
#' 
#' @return original SummarizedExperiment as in input
#' @return $output: list of pairwise partial correlation coefficients and pvalues, as well as the corresponding variable names
#' 
#' @examples
#' ... %>%
#'   mt_stats_multiv_net_GeneNet(name ="pcor") %>%
#' ...
#' 
#' @author EB
#' 

mt_stats_multiv_net_GeneNet = function(
  D,                       # SummarizedExperiment input
  name                     # unique name for this particular partial correlation matrix
) {
  
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))
  
  ## name
  if(missing(name))
    stop("name must be given")
  ## check for NA and throw an error if yes
  if(any(is.na(X)))
    stop("the data matrix contains NAs")
  
  require(GeneNet)
  require(igraph)
  
  # compute partial correlation using GeneNet
  pcor_GeneNet <- ggm.estimate.pcor(as.matrix(X), method = "dynamic")
  pval_GeneNet <- network.test.edges(pcor_GeneNet)
  
  # create result variables
  node1 <- colnames(pcor_GeneNet)[pval_GeneNet$node1]
  node2 <- colnames(pcor_GeneNet)[pval_GeneNet$node2]
  var <- paste0(node1,"_",node2, sep="")
  
  # create result table
  tab <- data.frame("var"=var, "statistic"=pval_GeneNet$pcor, "p.value"=pval_GeneNet$pval, "var1"=node1, "var2"=node2)
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = 'GeneNet partial correlation',
      output = list(
        table = tab,
        name = name,
        lstobj = NULL
        )
    )
  
  # return
  D
  
}