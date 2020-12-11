#' Scale data, mean 0 / sd 1 by default
#'
#' {ADD DESCRIPTION}
#'
#' @param D \code{SummarizedExperiment} input
#' @param center T/F, mean-center data? default: T
#' @param scale T/F, scale data to sd 1? default: T
#' @param ref_samples term which samples to use for center and scale calculcation
#'
#' @return assay: scaled data
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_trans_scale() %>% ...    # standard call, center and scale
#' ... %>% mt_pre_trans_scale(scale=F) %>% ...    # only mean centering
#' }
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_pre_trans_scale <- function(
  D,        # SummarizedExperiment input
  center=T, # mean 0?
  scale=T,   # SD 1?
  ref_samples #
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(is.logical(center))
  stopifnot(is.logical(scale))

  # scale
  if(!missing(ref_samples)){

    # get filtered samples
    filter_q <- dplyr::enquo(ref_samples)
    num_samp <- ncol(D)
    Da <- D %>% mti_format_se_samplewise()
    samples.used <- mti_filter_samples(Da, filter_q, num_samp)

    # get and filter
    Da_filtered <- t(assay(D))[samples.used,]

    if(center == T){
      # center by mean of selected samples
      da_f_means <- apply(Da_filtered, 2, mean, na.rm=T)
      Da <- sweep(Da, 2, da_f_means, FUN = "-")
    }
    if(scale == T){
      # scale by sd of selected samples
      da_f_sd <- apply(Da_filtered, 2, stats::sd, na.rm=T)
      Da <- sweep(Da, 2, da_f_sd, FUN = "/")
    }

    assay(D) = t(Da)

  } else {
    # if no sample filter given, just use scale() function
    assay(D) = t(scale(t(assay(D)),center=center,scale=scale))
  }

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('scaled, center=%d, scale=%d', center, scale),
      logshort = 'scaled'
    )

  # return
  D

}
