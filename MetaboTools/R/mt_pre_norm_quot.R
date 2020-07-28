#' Quotient normalization
#'
#' Implementation according to Dieterle et al., 2006\cr
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/16808434}{https://www.ncbi.nlm.nih.gov/pubmed/16808434}
#'
#' @param D \code{SummarizedExperiment} input
#' @param vars index vector of variables of be used, default: all
#' @param na_err T/F, throw error for NA's or just ignore?
#' @param ref_samples expression filtering reference samples from colData
#' @param met_max maximum fraction of missingness to select metabolites to be used in the reference (default: 1, i.e. all metabolites)
#'
#' @return assay: quotient-normalized version
#' @return $result: dilution factors
#'
#' @examples
#' \dontrun{#' # in the context of a SE pipeline
#' ... %>% mt_pre_norm_quot() %>% ...    # standard call
#' ... %>% mt_pre_norm_quot(ref_samples = GROUP=="ctrl") %>% ...    # use reference samples where 'GROUP' field in colData is 'ctrl'
#' ... %>% mt_pre_norm_quot(met_max = 0.2) %>% ...    # use only metabolites with <= 20% missing values to compute the reference used for normalization
#' }
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export
mt_pre_norm_quot <- function(
  D,                 #
  vars=1:dim(D)[1],  #
  na_err=F,         #
  ref_samples=NULL,   #
  met_max=1           #
) {

  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(!(met_max<0 || met_max>1))
  X = t(assay(D))
  if (any(unlist(X)[!is.na(unlist(X))]<0)) stop("Matrix contains negative values. Did you input logged data?")

  # check if there are any NAs
  if (sum(is.na(X[,vars]))>0) {
    # throw warning or error?
    if (na_err) {
      stop('Data matrix contains NAs')
    } else {
      MetaboTools:::mti_logwarning('Data matrix contains NAs')
    }
  }

  # apply reference sample filter?
  if (!missing(ref_samples)) {
    # select those samples according to the expression
    sample_filter_q <- dplyr::enquo(ref_samples)
    inds <- colData(D) %>%
      as.data.frame() %>%
      dplyr::mutate(tmporder=1:ncol(D)) %>%
      dplyr::filter(!!sample_filter_q) %>%
      .$tmporder
    useref <- rep(F, ncol(D))
    useref[inds] <- T
    # if no samples are left, throw error
    if (sum(useref) == 0) stop(sprintf("No samples match filter for quotient normalization: %s", dplyr::quo_name(sample_filter_q)))
  } else {
    useref = rep(T, ncol(D))
  }

  # compute metabolite missingness rate
  metmiss <- sapply(1:dim(X)[2], function(k){
    sum(is.na(X[,k]))/dim(X)[1]
  })
  # filter metabolites with missingness rate greater than met_max
  vars <- vars[metmiss <= met_max]
  # median reference sample
  ref = apply(X[useref,vars],2,function(x)stats::median(x,na.rm=T))
  # get dilution factors
  d = apply(X[,vars],1,  function(s) stats::median(as.numeric(s/ref),na.rm=T))
  # apply to each sample  (for each row=sample, divide values by median dilution factor)
  Y = t(sapply(1:dim(X)[1], function(i)unlist(X[i,]/d[i])))

  #Y = t(apply(X,1,  function(s) s /  d) )
  rownames(Y) = rownames(X)

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue::glue('quotient normalization based on {sum(useref)} reference samples and {length(vars)} variables: {dplyr::enquo(ref_samples) %>% as.character()}'),
      output = list(dilution=d)
    )

  # return
  assay(D) = t(Y)
  D

}
