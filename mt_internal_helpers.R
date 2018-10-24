# # MetaboTools
#
# Helper functions.
#
# last update: 2018-10-11
# authors: JK
#

# safely add to end of list, even if list does not exist
mti_add_to_list <- function(lst, element) {
  # init if needed
  if (is.null(lst) || length(lst)==0) lst=c()
  # add and return
  lst[[length(lst)+1]] = element
  lst
}

# either throws an error if an annotation field does not exist, or returns the annotation vector
# can be set to throw error if argument is not factor or vector of strings
mti_get_sample_anno <- function(D, strarg, requireFactor=F) {
  if (!(strarg %in% colnames(colData(D)))) stop(sprintf("'%s' not found in sample annotations",strarg))
  v = colData(D)[[strarg]]
  if (requireFactor && !(is.factor(v) || is.character(v))) stop(sprintf("'%s' is neither factor or vector of strings",strarg))
  v
}

# returns data frame as samples X variables, merges all sample annotations, and adds sample rownames as "primary" field
mti_format_se_samplewise <- function(se){
  cbind(colData(se),
        t(assay(se))) %>%
    as.data.frame() %>%
    rownames_to_column("primary")
}


# helper to generate metadata(.)$results entry
# will automatically generate fun and args fields, and send logtxt through logmsg()
mti_generate_result <- function(
  lst,              # output list, function should be used like this:  metadata(D)$results %<>% mti_generate_result(...)
  funargs,          # output from mti_funargs() that should be collected by calling function
  logtxt="",        # log text describing what the function did
  logshort=logtxt,  # shorter version of log text for printing; especially important for preprocessing funs... if not given, logtxt will be used
  output=NULL       # output structure of that function (e.g. plot, statistical results)... default is NULL (i.e. no output)
  ) {
  
  # ensure structure of funargs
  stopifnot("fun" %in% names(funargs))
  stopifnot("args" %in% names(funargs))
  
  # assemble list
  mti_add_to_list(
    lst,
    list(
      fun=funargs$fun,
      args=funargs$args,
      logtxt=logmsg(logtxt),
      logshort=logshort,
      output=output
    )
  )
}

# returns list of all arguments supplied to function 
# removes first argument $D if it exists
# remove mt from function name list
# returns list of $fun, already exploded, and $args
mti_funargs <- function(...) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))
  
  for(i in setdiff(names(formals), names(call)))
    call[i] <- list( formals[[i]] )
  
  # assemble results
  raw <-  as.list(match.call(sys.function(sys.parent(1)), call))
  res <- list(
    fun = strsplit(as.character(raw[[1]]), '_')[[1]],
    args = raw[-1]
  )
  # remove "mt"
  res$fun = res$fun[res$fun!="mt"]
  # make sure D does not exist in args
  res$args$D = NULL
  # return 
  res
 
}


# turns ... argument into key1=value1, key2=value2 etc.
mti_dots_to_str <- function(...) {
  l = eval(substitute(alist(...)))
  paste(sapply(names(l), function(k){sprintf('%s=%s',k,as.character(l[[k]]))}), collapse = ', ')
}

# generate extractor for $results entry, filter on top level function name
mti_res_get_toplevel <- function(D,str){metadata(D)$results[sapply(sapply(metadata(D)$results,'[[', 'fun',simplify=F), '[[', 1,simplify=F)==str]}
# extracts all preprocessing $results entries
mti_res_get_pre_entries  <- function(D){mti_res_get_toplevel(D,"pre")}
# extracts all plotting $results entries
mti_res_get_plot_entries <- function(D){mti_res_get_toplevel(D,"plots")}
# extracts all plotting $results entries
mti_res_get_plot_entries <- function(D){mti_res_get_toplevel(D,"stats")}
# extract plots only, not entire $results structure
mti_res_get_plots <- function(D,unlist=T){
  l=sapply(mti_res_get_toplevel(D,"plots"),'[[','output',simplify=F)
  if(unlist) l <- unlist(l,recursive=F)
  l
}

