# # MetaboTools
#
# Helper functions.
#
# last update: 2018-10-11
# authors: JK
#

# safely add to end of list, even if list does not exist
add_to_list <- function(lst, element) {
  # init if needed
  if (is.null(lst) || length(lst)==0) lst=c()
  # add and return
  lst[[length(lst)+1]] = element
  lst
}

# either throws an error if an annotation field does not exist, or returns the annotation vector
# can be set to throw error if argument is not factor or vector of strings
get_sample_anno <- function(D, strarg, requireFactor=F) {
  if (!(strarg %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations",strarg))
  v = colData(D)[[strarg]]
  if (requireFactor && !(is.factor(v) || is.character(v))) stop(sprintf("'%s' is neither factor or vector of strings",strarg))
  v
}

# returns data frame as samples X variables, merges all sample annotations, and adds sample rownames as "primary" field
format_se_samplewise <- function(se){
  cbind(colData(se),
        t(assay(se))) %>%
    as.data.frame() %>%
    rownames_to_column("primary")
}


# like match.call(), but also returns default values
match.call.defaults <- function(...) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))
  
  for(i in setdiff(names(formals), names(call)))
    call[i] <- list( formals[[i]] )
  
  
  match.call(sys.function(sys.parent()), call)
}