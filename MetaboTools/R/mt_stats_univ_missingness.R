#' Perform missingness significance analysis.
#'
#' This function will determine if NAs significantly accumulate in one of the sample groups. It is recommended that this function is run without prior missing value filtering.
#'
#' @param D \code{SummarizedExperiment} input
#' @param comp_col sample annotation (colData) column to compare against
#' @param stat_name name of comparison for later reference
#' @param sample_filter sample filter term to restrict to.
#'
#' @return $result: statistics object
#'
#' @examples
#' \dontrun{# run on sample field 'Group', name output stats object 'miss'
#' ... %>% mt_stats_univ_missingness(comp_col = 'Group', stat_name='miss') %>% ...
#' }
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_stats_univ_missingness <- function(
  D,      # SummarizedExperiment input
  comp_col,   # sample annotation (colData) column to compare against
  stat_name,    # name of comparison
  sample_filter
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(length(comp_col)==1)

  # merge data with sample info
  Ds <- D %>% MetaboTools:::mti_format_se_samplewise()

  ## FILTER SAMPLES
  if(!missing(sample_filter)) {

    filter_q <- dplyr::enquo(sample_filter)
    num_samp <- ncol(Ds)
    samples.used <- MetaboTools:::mti_filter_samples(Ds, filter_q, num_samp)
    Ds <- Ds[samples.used,]

  } else {
    samples.used = rep(T, ncol(Ds))
  }

  # get variable to compare to
  if (!(comp_col %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", comp_col))
  fixorder = function(x){o= unique(as.character(x)); gdata::reorder.factor(x, new.order=o)}
  vc = fixorder(as.factor(colData(D)[[comp_col]]))
  if (length(levels(vc))<2) stop(sprintf("'%s' has less than 2 factor levels",comp_col))




  # run models
  rawres <- sapply(1:nrow(Ds), function(i){
    # for (i in 1:nrow(Ds)) {

    # get metabolite
    m <- assay(D)[i,]

    # construct table
    tab <- table(is.na(m), vc)
    # fix table
    if (!("TRUE" %in% rownames(tab))) {
      tab %<>% rbind(rep(0,ncol(tab)))
      rownames(tab)[2] = "TRUE"
    }
    if (!("FALSE" %in% rownames(tab))) {
      tab %<>% rbind(rep(0,ncol(tab)))
      rownames(tab)[2] = "FALSE"
    }
    # run test
    test <- stats::fisher.test(tab)
    if (!("estimate" %in% names(test))) test$estimate=NA
    # construct extra fields
    ex <- gdata::unmatrix(tab)
    names(ex) = gsub("TRUE","missing",gsub("FALSE","present", names(ex)))
    # return
    as.data.frame(cbind(data.frame(var=rownames(D)[i], statistic=test$estimate, p.value=test$p.value),t(as.data.frame(ex))))

    # }
  },simplify=F)
  # combine
  res <- do.call(rbind, rawres)
  rownames(res) <- NULL

  ## add status information & results
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("missingness analysis with variable %s", as.character(comp_col)),
      output = list(
        table   = res,
        name    = stat_name,
        samples.used = samples.used,
        outcome = comp_col
      )
    )

  ## return
  D

}



