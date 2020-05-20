library(gdata)

#' Perform missingness significance analysis.
#' 
#' This function will determine if NAs significantly accumulate in one of the sample groups. It is recommended that this function is run without prior missing value filtering.
#'
#' @param D \code{SummarizedExperiment} input
#' @param comp sample annotation (colData) column to compare against
#' @param name name of comparison for later reference
#' @param samplefilter sample filter term to restrict to. 
#'
#' @return $result: statistics object
#' 
#' @examples
#' # run on sample field 'Group', name output stats object 'miss'
#' ... %>% mt_stats_univ_missingness(comp = 'Group', name='miss') %>% ...
#' 
#' @author JK
#' 
mt_stats_univ_missingness <- function(
  D,      # SummarizedExperiment input
  comp,   # sample annotation (colData) column to compare against
  name,    # name of comparison
  samplefilter
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(length(comp)==1)
  
  # filter samples?
  Ds <- D
  if(!missing(samplefilter)) {
    # translate lazy expression
    filter_q <- enquo(samplefilter)
    # get indices
    keep <-     
      colData(D) %>% as.data.frame() %>%
      add_rownames() %>%
      filter(!!filter_q) %>%
      `[[`("rowname") %>%
      as.numeric()
    # filter SE
    Ds <- Ds[,keep]
  }
  
  # get variable to compare to
  if (!(comp %in% colnames(colData(Ds)))) stop(sprintf("'%s' not found in sample annotations.", comp))
  fixorder = function(x){o= unique(as.character(x)); gdata::reorder.factor(x, new.order=o)} 
  vc = fixorder(as.factor(colData(Ds)[[comp]]))
  if (length(levels(vc))<2) stop(sprintf("'%s' has less than 2 factor levels",comp))
  
  
  
  
  # run models
  rawres <- sapply(1:nrow(Ds), function(i){
    # for (i in 1:nrow(Ds)) {
    
    # get metabolite
    m <- assay(Ds)[i,]
    
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
    test <- fisher.test(tab)
    if (!("estimate" %in% names(test))) test$estimate=NA
    # construct extra fields
    ex <- unmatrix(tab)
    names(ex) = gsub("TRUE","missing",gsub("FALSE","present", names(ex)))
    # return
    as.data.frame(cbind(data.frame(var=rownames(Ds)[i], statistic=test$estimate, p.value=test$p.value),t(as.data.frame(ex))))
    
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
      logtxt = sprintf("missingness analysis with variable %s", as.character(comp)),
      output = list(
        table   = res,
        name    = name
      )
    )
  
  ## return
  D
  
  
}



