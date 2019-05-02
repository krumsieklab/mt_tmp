# MetaboTools
#
# Univariate missingness analysis with a factor outcome.
#
# last update: 2019-01-06
# authors: JK
#

library(gdata)

mt_stats_univ_missingness <- function(
  D,      # SummarizedExperiment input
  comp,   # sample annotation (colData) column to compare against
  name,    # name of comparison
  samplefilter
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  stopifnot(length(comp)==1)
  
  # get variable to compare to
  if (!(comp %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations.", comp))
  vc = mti_fixorder(as.factor(colData(D)[[comp]]))
  if (length(levels(vc))<2) stop(sprintf("'%s' has less than 2 factor levels",comp))
  
  # filter samples?
  Ds <- D
  if(!missing(samplefilter)) {
    filter_q <- enquo(samplefilter)
    cd <- colData(D) %>% as.data.frame() %>% filter(!!filter_q)
    # keep <-
  }
  
  
  
  # run models
  rawres <- sapply(1:nrow(D), function(i){
    # for (i in 1:nrow(D)) {
    
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
    test <- fisher.test(tab)
    if (!("estimate" %in% names(test))) test$estimate=NA
    # construct extra fields
    ex <- unmatrix(tab)
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
      logtxt = sprintf("missingness analysis with variable %s", as.character(comp)),
      output = list(
        table   = res,
        name    = name
      )
    )
  
  ## return
  D
  
  
}



