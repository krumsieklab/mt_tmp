#' Computes partial correlation matrix using the GeneNet estimator.
#'
#' Implementation according to Schäfer and Strimmer, 2006\cr
#' https://www.ncbi.nlm.nih.gov/pubmed/16646851
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name name of the correlation matrix
#' @param sample_filter term defining which samples to use for GGM calculation (default: all samples)
#'
#' @return original SummarizedExperiment as in input
#' @return $output: list of pairwise partial correlation coefficients and pvalues, as well as the corresponding variable names
#'
#' @examples
#' \dontrun{... %>%
#'   mt_stats_multiv_net_GeneNet(stat_name ="pcor") %>%
#' ...}
#'
#' @author EB
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_stats_multiv_net_GeneNet = function(
  D,                       # SummarizedExperiment input
  stat_name,                    # unique name for this particular partial correlation matrix
  sample_filter
) {

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))

  ## stat_name
  if(missing(stat_name))
    stop("stat_name must be given")
  ## check for NA and throw an error if yes
  if(any(is.na(X)))
    stop("the data matrix contains NAs")

  ## FILTER SAMPLES?
  if(!missing(sample_filter)) {
    # merge data with sample info
    Ds <- D %>% mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK
    filter_q <- dplyr::enquo(samplefilter)
    Ds <- Ds %>%
      dplyr::mutate(tmpsamplenum = 1:nrow(Ds)) %>%
      dplyr::filter(!!filter_q) %>%
      droplevels()
    # message("filter metabolites: ", metab_filter_q, " [", nrow(stat), " remaining]")
    # did we leave 0 rows?
    if (nrow(Ds)==0) stop("Filtering left 0 rows")
    if (nrow(Ds)==ncol(D)) mti_logwarning('filtering did not remove any samples')

    # store used samples
    samples.used <- rep(F, ncol(D))
    samples.used[Ds$tmpsamplenum] <- T

  } else {
    samples.used = rep(T, ncol(D))
  }

  # filter
  X <- X[samples.used,]

  # compute partial correlation using GeneNet
  pcor_GeneNet <- GeneNet::ggm.estimate.pcor(as.matrix(X), method = "dynamic", verbose=FALSE)
  pval_GeneNet <- GeneNet::network.test.edges(pcor_GeneNet, plot=FALSE)

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
        name = stat_name,
        lstobj = NULL,
        outcome = NA
        )
    )

  # return
  D

}
