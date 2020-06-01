#' Filter by missingness.
#'
#' Filters either samples or metabolites, won't do both in one call.
#'
#' @param D \code{SummarizedExperiment} input
#' @param met_max Maximum fraction of missing metabolites (between 0 and 1.0)
#' @param sample_max Maximum fraction of missing samples (between 0 and 1.0)
#' @param met_group Name of of a colData sample annotation column; met_max will be applied to each group separately, metabolite must have at most met_max in any of the groups.
#' @param log_filtered Write list of filtered metabolites into the log text? Default: False
#' @param log_sample_col Required if log_filtered=T for sample filtering. Specifies which column of colData() to output into the log.
#'
#' @return SE rows or columns will be filtered.
#' @return $output: logical vector of kept metabolites/samples
#'
#' @examples
#' \dontrun{# first remove samples with >10% missingness, then metabolites with >20% missingness
#' ... %>%
#'   mt_pre_filtermiss(sample_max=0.1) %>%
#'   mt_pre_filtermiss(met_max=0.2) %>%
#' ...}
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_pre_filtermiss <- function(
  D,                   # SummarizedExperiment input
  met_max=NA,           # maximum fraction of missing metabolites
  sample_max=NA,        # " samples
  met_group = NA,      # for each group of samples met_max be applied
  log_filtered = F,    # write filtered ones into log?
  log_sample_col=''  # which field for sample anno?
) {

  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(!is.na(met_max) || !is.na(sample_max))
  stopifnot(!(is.na(met_max) && is.na(sample_max)))
  stopifnot(!(!is.na(met_max) && (met_max<0 || met_max>1)))
  stopifnot(!(!is.na(sample_max) && (sample_max<0 || sample_max>1)))

  # perform filtering
  if (!is.na(met_max)) {
    NA.mat = is.na(assay(D))
    # metabolite
    lst_filtered = c()
    if(is.na(met_group)){ # if group met_max
      na.stat = rowSums(NA.mat)
      metKeep = na.stat/ncol(D) <= met_max
      # store list of filtered
      lst_filtered = D %>% rowData() %>% .$name %>% .[!metKeep]
      # filter
      D=D[metKeep,]
      na.stat[metKeep]
    }else{
      stopifnot(met_group %in% (colData(D) %>% colnames))
      xmet_group = colData(D)[,met_group]
      na.stat = xmet_group %>% unique %>% {. = as.character(.); names(.) = .; .} %>%
        sapply(function(x) rowSums(NA.mat[, xmet_group == x])/sum(xmet_group == x))
      metKeep = rowSums( na.stat <= met_max ) == ncol(na.stat)
      # store list of filtered
      lst_filtered = D %>% rowData() %>% .$name %>% .[!metKeep]
      # filter
      D=D[metKeep,]
      na.stat[metKeep, ]
    }

    # generate logtxt
    logtxt = sprintf('metabolites filtered, %.2f%%, %d of %d removed', round(met_max*100),sum(!metKeep),length(metKeep))
    if (log_filtered) {
      logtxt %<>% paste0(sprintf("  -  list of filtered metabolites: %s", paste0(lst_filtered, collapse=" / ")))
    }

    # add status information
    funargs <- mti_funargs()
    metadata(D)$results %<>%
      mti_generate_result(
        funargs = funargs,
        logtxt = logtxt,
        logshort = sprintf("filter met %.2f%%", round(met_max*100)),
        output = list(kept=as.vector(metKeep), na.stat = na.stat, na.mat = NA.mat[metKeep,])
      )

  } else {

    # sample
    sampleKeep = apply(is.na(assay(D)),2,sum)/nrow(D) <= sample_max

    # generate logtxt
    logtxt = sprintf('samples filtered, %.2f%%, %d of %d removed', sample_max*100,sum(!sampleKeep),length(sampleKeep))
    if (log_filtered) {
      if (nchar(log_sample_col)==0) stop("If log_filtered=T and filtering samples, log_sample_col must be specified.")
      lst_filtered = D %>% colData() %>% .[[log_sample_col]] %>% .[!sampleKeep]
      logtxt %<>% paste0(sprintf("  -  list of filtered metabolites: %s", paste0(lst_filtered, collapse=" / ")))
    }

    # filter
    D=D[,sampleKeep]

    # add status information
    funargs <- mti_funargs()
    metadata(D)$results %<>%
      mti_generate_result(
        funargs = funargs,
        logtxt = logtxt,
        logshort = sprintf("filter samples %.2f%%", sample_max*100),
        output = list(kept=as.vector(sampleKeep))
      )

  }

  # throw error if filtering caused empty dataset
  if (ncol(D)==0 || nrow(D)==0) stop("Filtering resulted in empty dataset.")

  # return
  D
}


