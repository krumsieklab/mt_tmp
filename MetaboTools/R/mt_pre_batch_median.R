#' Median batch correction.
#'
#' Same approach that Metabolon uses for runday correction. Divides metabolite values by the median value per batch.
#'
#' @param D \code{SummarizedExperiment} input
#' @param batches sample annotation (colData) column name that contains batch assignment
#' @param ref_samples expression to filter out reference samples to use from colData
#'
#' @return assay: batch-corrected version
#'
#' @examples
#' \dontrun{... %>% mt_pre_batch_median(batches="BATCH") %>% ...
#' }
#'
#' @author JK
#'
#' @export
mt_pre_batch_median = function(
  D,          # SummarizedExperiment input
  batches,    # sample annotation column that contains batch info
  ref_samples # expression to filter out reference samples to use from colData
) {

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  X = t(assay(D))

  # get variable
  if (!(batches %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", batches))
  b = colData(D)[[batches]]
  # ensure it's a factor or character vector
  if (!(is.character(b) || is.factor(b))) stop(sprintf("'%s' has to be character, factor, or numeric.", batches))
  b = as.factor(b)

  # no negative values allowed
  if (min(X,na.rm=T)<0) stop("Matrix contains negative values.")
  # check if data actually have been logged by preprocessing
  if (length(MetaboTools:::mti_res_get_path(D, c("pre","trans","log"))) > 0 |
      length(MetaboTools:::mti_res_get_path(D, c("flag", "logged"))) > 0)
    stop("Median batch correction can only be performed on non-logged data.")

  # get samples to use for median calculation
  if (!missing(ref_samples)) {
    # select those samples according to the expression
    ref_samples_q <- dplyr::enquo(ref_samples)
    inds <- colData(D) %>%
      as.data.frame() %>%
      dplyr::mutate(tmporder=1:ncol(D)) %>%
      dplyr::filter(!!ref_samples_q) %>%
      .$tmporder
    use_samples <- rep(F, ncol(D))
    use_samples[inds] <- T
  } else {
    # all samples
    use_samples <- rep(T, ncol(D))
  }

  # median per metabolite
  for (i in 1:length(levels(b))) {
    batch = levels(b)[i]
    # check that there are any reference samples available
    if (sum(b==batch & use_samples) == 0) stop(sprintf("No reference samples for batch '%s'", batch))
    # build median vector for this batch
    med <- X[b==batch & use_samples,,drop=F] %>% apply(2, stats::median, na.rm=T)
    # transform into matrix of the same size as the batch
    med_matrix <- replicate(sum(b==batch), med) %>% t()
    # median normalize
    X[b==batch,] <-  X[b==batch,] / med_matrix
  }

  # replace original assay with batch corrected assay
  assay(D) = t(X)

  # ref samples logging string
  refadd <- if(missing(ref_samples)){""}else{sprintf(": %s", as.character(dplyr::enquo(ref_samples)))}

  # add status information
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("median-scaling per batch in '%s', based on %d reference samples%s",batches,sum(use_samples),refadd)
    )

  # return
  D


}
