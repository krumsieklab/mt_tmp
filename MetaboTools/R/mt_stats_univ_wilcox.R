#' Wilcox test
#'
#'
#' @param D \code{SummarizedExperiment} input
#' @param y name of the column in colData to compare with metabolite
#' @param stat_name name under which this comparison will be stored, must be unique to all other statistical results
#' @param sample_filter term which samples to filter to first... e.g. used if the data contains >2 groups but the user wants to run a two-group comparison
#' @param exact_flag optional to set the exact flag in wilcox.test function
#' @param paired_flag optional to set the paired flag in wilcox.test function, would be applied for numeric y
#' @return $result: statistics object
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#' @import glue
#'
#' @examples
#' \donttest{# run lm with no confounders, "Group" as outcome
#' # filter to groups "Li_2" and "Li_5"
#' # name the comparison "Li's"
#' ... %>%
#'  mt_stats_univ_wilcox(
#'    y     = Group,
#'    sample_filter = (Group %in% c("Li_2","Li_5")),
#'    stat_name         = "Li's"
#'  ) %>% ...
#'  }
#'
#' @author RB
#'
#' @export

mt_stats_univ_wilcox <- function(
  D,              # SummarizedExperiment input
  y,        # name of the column in colData to compare with metabolite
  stat_name,      # name of comparison,
  sample_filter,
  paired_flag=F, # if tests should be paired
  exact_flag=NULL
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  # check that y is in the colData
  if (!(y %in% colnames(colData(D)))) stop(sprintf("There is no column called %s in the colData", y))
  # make sure name does not exist yet
  if (stat_name %in% unlist(MetaboTools:::mti_res_get_stats_entries(D) %>% purrr::map("output") %>% purrr::map("name"))) stop(sprintf("stat element with name '%s' already exists",stat_name))
  
  # merge data with sample info
  Ds <- D %>% MetaboTools:::mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 8/17/20, JK
  
  ## FILTER SAMPLES
  if(!missing(sample_filter)) {
    
    filter_q <- dplyr::enquo(sample_filter)
    Ds <- Ds %>%
      dplyr::mutate(tmpsamplenum = 1:nrow(Ds)) %>%
      dplyr::filter(!!filter_q) %>%
      droplevels()
    # message("filter metabolites: ", metab_filter_q, " [", nrow(stat), " remaining]")
    # did we leave 0 rows?
    if (nrow(Ds)==0) stop("Filtering left 0 rows")
    if (nrow(Ds)==ncol(D)) mti_logwarning('filtering did not filter out any samples')
    
    # store used samples
    samples.used <- rep(F, ncol(D))
    samples.used[Ds$tmpsamplenum] <- T
    # drop dummy column
    Ds %<>% dplyr::select(-tmpsamplenum)
    
  } else {
    samples.used = rep(T, ncol(D))
  }
  
  # check that outcome is either binary or numerical and run test accordingly
  mets <- rownames(D)
  outvec <- Ds[[y]]
  cl <- outvec %>% class()
  
  if (("character" %in% cl) || ("factor" %in% cl)) {
    if ((outvec %>% as.factor() %>% levels() %>% length()) != 2) {
      stop("If outcome is a factor, it must have exactly two levels")
    }
    # save groups
    outgroups <- outvec %>% as.factor() %>% levels()
    # run wilcox
    # we want to model metabolite ~ outcome
    wt <- lapply(mets, function(x){
      input <- split(Ds[,x], Ds[[y]])
      res  <- wilcox.test(input[[1]], input[[2]], alternative = "two.sided", exact=exact_flag)
      res <- data.frame("statistic"=res$statistic, "p.value"=res$p.value, "method"=res$method)
      return(res)
    }) %>% do.call(rbind,.) %>% data.frame()

  } else {
    outgroups <- NULL
    # run wilcox
    # we want to model metabolite ~ outcome
    wt <- lapply(mets, function(x){
      res  <- wilcox.test(Ds[,x], Ds[[y]], paired=paired_flag, alternative = "two.sided", exact=exact_flag)
      res <- data.frame("statistic"=res$statistic, "p.value"=res$p.value, "method"=res$method)
      return(res)
    }) %>% do.call(rbind,.) %>% data.frame()
  }
  
  # add columns with metabolite names and y variable 
  wt %<>% mutate(var=mets, term=y)
  # reorder columns
  wt %<>% select(var, term, statistic, p.value)
  # rearrange back to original metabolite order
  o <- match(rownames(D),wt$var)
  stopifnot(!any(is.na(o))) # sanity check
  wt <- wt[o,]
  # make sure that NAs in the outcome are set to FALSE in the list of used samples
  samples.used[is.na(Ds[[y]])] <- F
  
  ## add status information & results
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Wilcox rank sum test, %s", y),
      output = list(
        table = wt,
        #formula = as.formula(glue::glue("~ {y}")),
        name    = stat_name,
        groups = outgroups,
        #lstobj = NULL,
        samples.used = samples.used,
        outcome = y
      )
    )
  
  ## return
  D
  
}
