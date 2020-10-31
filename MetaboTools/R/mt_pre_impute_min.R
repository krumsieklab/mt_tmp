#' Impute using minimum value of that metabolite across the samples
#'
#'
#' @param D \code{SummarizedExperiment} input
#' @param verbose T/F, whether to output intermediate steps, default: F
#' @return assay: imputed data
#'
#' @examples
#' \dontrun{# in the context of a SE pipeline
#' ... %>% mt_pre_impute_min() %>%
#' }
#'
#' @author RB
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_pre_impute_min <- function(
  D,      # SummarizedExperiment input
  verbose=F
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  my_min <- function(x) {
    if(all(is.na(x))){
      return(NA)
    } else{
      return(min(x, na.rm=T))
    }
  }
  
  df <- D %>% assay() %>%  data.frame()
  incom.obs <- which(apply(df,2,function(x) any(is.na(x))))
  if(verbose)message(paste0("Number of imcomplete observations: ", length(incom.obs)))
  all_nas <- length(which(apply(df,1,function(x) all(is.na(x)))))
  # impute
  df = apply(df, 1, function(x) {x[is.na(x)] <-  my_min(x); x} ) %>% t()
  assay(D) = df

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('imputed via minimum value, %d metabolites with all NAs, returned as NAs', all_nas),
      logshort = 'impute min'
    )

  # return
  D

}
