# MetaboTools
#
# Aggregate pathway values. Takes a pathway annotation column of the metabolites 
# and builds one representative sample per pathway for each sample.
#
# Implemented approaches:
# 1. Eigenmetabolite/eigengene/eigenvalue PCA based approach. Data matrix cannot have NAs.
# 2. Mean scaled value. Will scale() data before. Data matrix can have NAs.
#
# last update: 2018-10-11
# authors: JK
#

mt_modify_aggPW <- function(
  D,          # SummarizedExperiment input
  pw,         # string, column name for pathway annotations
  method      # one of: "eigen", "aggmean"
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  if (!(method %in% c("eigen","aggmean"))) stop("'method' must be either 'eigen' or 'aggmean'")
  
  # get variable
  if (!(pw %in% colnames(rowData(D)))) stop(sprintf("'%s' not found in metabolite annotations.", pw))
  p = rowData(D)[[pw]]
  # ensure it's a factor or character vector
  if (!(is.character(p) || is.factor(p))) stop(sprintf("'%s' has to be character, factor, or numeric.", pw))
  p = mti_fixorder(as.factor(p))
  
  # agg calculation
  X = t(assay(D))
  output <- NULL
  if (method=="eigen") {
    # for eigenvalue methods, matrix cannot have any NAs
    if (any(is.na(X))) stop("no NA values allowed for 'eigen' method")
    # calc
    res <- unique(p) %>% mti_removeNAs() %>% lapply(function(g){
      met <- mti_removeNAs(p==g, replaceWith=F)
      pca = prcomp(as.data.frame(X[,met]))
      list(pc1=pca$x[,1],
           expvar=(pca$sdev)^2 / sum(pca$sdev^2) )
    })
    # assemble into matrix
    M = as.data.frame(sapply(res, function(x){x$pc1}))
    colnames(M) <- unique(p) %>% mti_removeNAs()
    # return explained variances in output
    output$expvar = map(res, "expvar")
    
  } else if (method=="aggmean") {
    # TODO implement me
    # result has to go into M
    
    M <- unique(p) %>% mti_removeNAs() %>% sapply(function(g){
      met <- mti_removeNAs(p==g, replaceWith=F)
      apply(as.data.frame(X[,met]),1,mean, na.rm=T  )
    })
    colnames(M) <- unique(p) %>% mti_removeNAs()
    
  } else {stop("bug")}
  
  # prepare matrix
  full.names <- colnames(M)
  colnames(M) <- make.names(full.names)
  
  # check which variables can be copied [all of this can probable be done simpler]
  copyover <- sapply(colnames(rowData(D)), function(c) {
    # verify variable, there must be only one value for each instance
    all( (rowData(D) %>% as.data.frame() %>% as.tibble() %>% group_by_(pw) %>% summarise(n_distinct(!!rlang::sym(c))))[[2]] ==1 )
  })
  # generate new rowData
  rd <- rowData(D) %>% as.data.frame() %>% as.tibble() %>% filter(!duplicated(!!rlang::sym(pw))) %>% 
    select(which(copyover)) %>% filter(!is.na(!!rlang::sym(pw))) %>% mutate(name=!!rlang::sym(pw)) %>% select(name,everything())
  # sanity bug check
  stopifnot(all.equal(rd$name, rd[[pw]]))
  
  # generate summarized experiment
  newD <- SummarizedExperiment(assay = t(M),
                               colData  = colData(D),
                               rowData  = rd,
                               metadata = metadata(D)
  )
  
  # add status information
  funargs <- mti_funargs()
  metadata(newD)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("pathway aggregation by '%s', %d pathway scores generated", pw, nrow(newD))
    )
  
  # return
  newD
  
  
}
