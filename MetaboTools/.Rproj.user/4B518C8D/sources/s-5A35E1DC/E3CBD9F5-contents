library(stats)

#' Computes Kendall's rank correlation.
#' If present, NAs will be omitted.
#' 
#'
#' @param D \code{SummarizedExperiment} input
#' @param var string name of the colData variable to use as ordered categorical variable for the correlation calculation. class(D[[var]]) needs to be numeric.
#' @param name name under which this comparison will be stored, must be unique to all other statistical results
#' @param samplefilter optional sample filter condition
#' @param exactFlag optional to set the exact flag in cor.test function
#' 
#' @return original SummarizedExperiment as in input
#' @return $output: list of Kendall's correlation coefficients and pvalues, as well as the corresponding variable names
#' 
#' @examples
#' ... %>%
#'   mt_stats_univ_tau(var = "Stage", samplefilter = (GROUP %in% "Tumor"), name = "tau") %>%
#' ...
#' 
#' @author EB, RB (modified on 2020-05-21)
#' 
#' @export

mt_stats_univ_tau = function(
  D,
  var,
  name,
  samplefilter, exactFlag=NULL){
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  # check that "var" is in the colData
  if (!(var %in% colnames(colData(D)))) stop(sprintf("There is no column called %s in the colData", var))
  # check that "var" is numeric
  if (class(D[[var]])!="numeric") stop(sprintf("%s must be numeric",var))
  
  # make sure name does not exist yet
  if (name %in% unlist(mti_res_get_stats_entries(D) %>% map("output") %>% map("name"))) stop(sprintf("stat element with name '%s' already exists",name))
  
  # merge data with sample info
  Ds <- D %>% mti_format_se_samplewise() 
  
  ## FILTER SAMPLES
  if(!missing(samplefilter)) {
    
    filter_q <- enquo(samplefilter)
    Ds <- Ds %>%
      filter(!!filter_q) %>%
      droplevels()
    # message("filter metabolites: ", metab_filter_q, " [", nrow(stat), " remaining]")
    # did we leave 0 rows?
    if (nrow(Ds)==0) stop("Filtering left 0 rows")
    if (nrow(Ds)==ncol(D)) mti_logwarning('filtering did not filter out any samples')
    
  }

  met <- colnames(Ds)[(length(colnames(colData(D)))+2):length(colnames(Ds))]
  # compute association to the phenotype
  rr <- lapply(met, function(x){
    d=cor.test(Ds[,x], Ds[[var]], method="kendall", alternative = "two.sided")
    list("statistic"=d$estimate, "p.value"=d$p.value, exact=exactFlag)
  })
  names(rr) <- met
  
  # revert list structure
  revert_list_str <- function(ls) {
    # get sub-elements in same order
    x <- lapply(ls, `[`, names(ls[[1]]))
    # stack and reslice
    apply(do.call(rbind, x), 2, as.list) 
  }
  rr_reverse <- revert_list_str(rr)
  
  # arrange results in dataframe
  tab <-cbind.data.frame(as.data.frame(do.call(rbind, rr_reverse$statistic)),
                                   as.data.frame(do.call(rbind, rr_reverse$p.value)))
  colnames(tab) <- c("statistic","p.value")
  # add term column with ordinal variable
  tab$term <- rep(var, dim(tab)[1])
  # add column with names
  tab$var <- rownames(tab)
  
  ## construct output groups variable
  outgroups <- unique(Ds[[var]])
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = 'Kendall rank correlation (tau)',
      output = list(
        table = tab,
        name = name,
        lstobj = NULL,
        groups = outgroups
      )
    )
  
  # return
  D
  
}