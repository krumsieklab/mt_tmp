# # MetaboTools
#
# Helper functions.
#
# last update: 2020-09-17
# authors: JK,MB, KC
#

#' Return stats output by name.
#'
#' Finds a named statistical result from a mt_stats... function and returns the $output$table dataframe.
#'
#' @param D SummarizedExperiment
#' @param name Name of statistical comparison
#' @param fullstruct optional, output entire $output structure, not just the $table inside.
#'
#' @return $output$table dataframe
#' @noRd
mti_get_stat_by_name <- function(D, name, fullstruct=F){
  stopifnot("SummarizedExperiment" %in% class(D))

  if(! ("results" %in% names(metadata(D))))
    stop("no results element found in D")

  stats <- mti_res_get_stats_entries(D)

  if(length(stats) == 0)
    stop("no stats element found in D")

  names  <- stats %>% purrr::map_chr(~.x$output$name)
  output <- which(names == name)

  if(length(output) == 0)
    stop("stat element with name ", name, " does not exist")
  if(length(output)  > 1)
    stop("there are multiple stat elements with name ", name)

  if(!fullstruct) {
    output <- stats[[ output ]]$output$table
    if( ! any(c("tibble", "data.frame") %in% class(output)) )
      stop("output of stat result ", stat, " is not a table")

    if( ! ("var" %in% colnames(output)) )
      stop("output of stat result ", name, " does not have 'var' column")
  } else {
    output <- stats[[ output ]]$output
  }
  output
}

#' Retrieve ML model results by name
#'
#' Returns a list of ML results given a result name
#'
#' @param D SummarizedExperiment object
#' @param name the name of the ML result of interest
#'
#' @result res_list: list containing the output and output2 lists for a given name
#'
#' @noRd
mti_get_ml_res_by_name <- function(D, name){
 # code is identical to mti_get_stat_by_name, but returns output and output2 lists
  stopifnot("SummarizedExperiment" %in% class(D))

  if(! ("results" %in% names(metadata(D)))){
    stop("no results element found in D")
  }

  stats <- MetaboTools:::mti_res_get_path(D,"ml")

  if(length(stats) == 0){
    stop("no stats element found in D")
  }

  names  <- stats %>% purrr::map_chr(~.x$output$name)
  output <- which(names == name)

  if(length(output) == 0){
    stop("stat element with name ", name, " does not exist")
  }
  if(length(output)  > 1){
    stop("there are multiple stat elements with name ", name)
  }

  res_list <- list()
  res_list$output <- stats[[ output ]]$output
  res_list$output2 <- stats[[ output ]]$output2

  res_list

}


#' Safely add to end of list, even if list does not exist
#'
#' Adds an item to the end of a list, even if the list is empty or NULL
#'
#' @param lst list to be extended
#' @param element element to be added
#' @param oname name of the new element
#'
#' @return extended list
#'
#' @noRd
mti_add_to_list <- function(lst, element, oname) {
  # init if needed
  if (is.null(lst) || length(lst)==0) lst=c()
  # add and return
  lst[[length(lst)+1]] = element
  names(lst)[length(lst)] = oname

  lst
}


#' Concatenate data matrix with sample annotations
#'
#' Returns data frame as samples X variables, merges all sample annotations, and adds sample rownames as "merge.primary" field
#'
#' @param D SummarizedExperiment input
#'
#' @return Concatenated dataframe
#'
#' @noRd
mti_format_se_samplewise <- function(D){
  # coldata and assay cannot have overlapping names
  # (this should be caugh earlier in the pipeline, but here is where it causes trouble)
  inters <- intersect(colnames(colData(D)), rownames(D))
  if (length(inters)>0) {
    stop(sprintf("There are metabolites and colData variables with the same name: %s", paste0(inters, collapse = ", ")))
  }
  # cbind
  cbind(colData(D),
        t(assay(D))) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("merge.primary")
}


#' Helper function to generate metadata(.)$results entry
#'
#' Automatically generates fun and args fields, and sends logtxt through logmsg()
#'
#' @param lst result list, function should be used like this:  metadata(D)$results %<>% mti_generate_result(...)
#' @param funargs output from mti_funargs() that should be collected by calling function
#' @param logtxt log text describing what the function did
#' @param logshort shorter version of log text for printing; especially important for preprocessing funs... if not given, logtxt will be used
#' @param output output structure of that function (e.g. plot, statistical results)... default is NULL (i.e. no output)
#' @param output2 optional second output... default is NULL
#'
#' @return list ready to be stored in metadata()$results
#'
#' @noRd
mti_generate_result <- function(
  lst,              #
  funargs,          #
  logtxt="",        #
  logshort=logtxt,  #
  output=NULL,      #
  output2=NULL     #
) {

  # ensure structure of funargs
  stopifnot("fun" %in% names(funargs))
  stopifnot("args" %in% names(funargs))

  this.uuid = uuid::UUIDgenerate()

  # assemble list
  mti_add_to_list(
    lst,
    list(
      fun=funargs$fun,
      args=funargs$args,
      logtxt=mti_logmsg(logtxt),
      logshort=logshort,
      uuid=this.uuid,
      output=output,
      output2=output2
    ),
    oname = paste(paste(funargs$fun,collapse = "_"), this.uuid, sep = ".")
  )
}


#' Return list of all arguments supplied to function
#'
#' Removes first argument $D if it exists.
#' Removes "mt" from function name list
#' Returns list of $fun, already exploded (e.g. c("plots","boxplot")), and $args
#'
#' @param ... Not actually used (TODO delete?). Accesses arguments of parent function
#'
#' @return List of arguments, see description.
#'
#' @noRd
mti_funargs <- function(...) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))

  for(i in dplyr::setdiff(names(formals), names(call)))
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


#' Turns ... argument into string, key1=value1, key2=value2 etc.
#'
#' @param ... Arbitary list of input arguments
#'
#' @return String key/value list.
#'
#' @noRd
mti_dots_to_str <- function(...) {
  l = eval(substitute(alist(...)))
  paste(sapply(names(l), function(k){sprintf('%s=%s',k,as.character(l[[k]]))}), collapse = ', ')
}


#' Confounder correction for a SummarizedExperiment
#'
#' Used specifically for boxplot function to generate confounder-corrected residuals for plotting.
#'
#' @param D \code{SummarizedExperiment} input
#' @param formula formula for correction
#'
#' @returns SummarizedExperiment with corrected data
#'
#' @noRd
mti_correctConfounder <- function(D, formula){
  d <- D %>% mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK
  d_cor <- rownames(D) %>%
    purrr::map_dfc(function(m){
      f   <- stats::update.formula(formula, stringr::str_c(m, "~."))
      mod <- stats::lm(f, data = d, na.action = na.exclude)
      res <- stats::resid(mod)
      res
    }) %>%
    stats::setNames(rownames(D)) %>%
    as.matrix() %>% t()
  colnames(d_cor) <- colnames(D)
  assay(D)        <- d_cor
  D
}


#' Extract results from a "path" of entries in metadata
#'
#' Extracts metadata()$results entries in a given namespace of arbitrary depth. See examples below.
#'
#' @param D \code{SummarizedExperiment} input
#' @param path path, list of strings
#'
#' @returns SummarizedExperiment with corrected data
#'
#' @examples
#' \dontrun{# get all result entries for plots
#' D %>% mti_res_get_path("plots")
#' # get all result entries for stats
#' D %>% mti_res_get_path("stats")
#' # check if quotient normalization has been performed
#' q <- D%>% mti_res_get_path(c("pre","norm","quot"))
#' if (length(q)==0) stop("No quotient normalization performed.")
#' }
#'
#' @author JK
#' @noRd
mti_res_get_path <- function(D,path) {
  # [sorry for the code, but it works]s
  labels <- lapply(metadata(D)$results,'[[', 'fun')
  m <- rep(T,length(labels))
  # exclude non-matches
  for (i in 1:length(path)) {
    m <- m & sapply(labels, function(label){
      if (length(label)<i) F
      else {
        if (label[i]!=path[i])F
        else T
      }
    })
  }
  # return
  metadata(D)$results[m]
}


#' Extract all stats entries.
#'
#' Returns all stats entries from the $results list.
#'
#' @param D \code{SummarizedExperiment} input
#'
#' @return list if stats results
#'
#' @author JK
#' @noRd
mti_res_get_stats_entries <- function(D){mti_res_get_path(D,"stats")}


#' Add left- and right-aligned x axis labels to ggplot
#'
#' ggplot is missing the functionality to add x-axis labels that are left- and right-aligned. This function adds those.
#'
#' @param ggplot object
#' @param left text on the left
#' @param right text on the right
#'
#' @examples
#' \dontrun{# To demonstrate function behavior on empty plot:
#' (ggplot(mapping = aes(x=0,y=0)) + geom_point()) %>% mti_add_leftright_gg("left label","right label")
#' }
#'
#' @return ggplot object
#'
#' @author MB, JK
#' @noRd
mti_add_leftright_gg <- function(gg, left, right) {

  # with backward compatibility
  if (utils::compareVersion(as.character(utils::packageVersion("ggplot2")),"3.3.0")>=0) { # at least 3.3.0
    # ggplot2 version >= 3.3.0
    ggbld <- ggplot2::ggplot_build(gg)
    xticks =  ggbld$layout$panel_params[[1]]$x$get_breaks() # needed for >=3.3.0
    xticks.minor = ggbld$layout$panel_params[[1]]$x$get_breaks_minor() # needed for >=3.3.0
    xlims = ggbld$layout$panel_params[[1]]$x.range

    # add positions of annotation labels
    xticks = c(xlims[1], xticks, xlims[2])
    # get breaks labels
    xtlabs = ggbld$layout$panel_params[[1]]$x$get_labels() # needed for >=3.3.0

  } else {
    # ggplot2 version < 3.3.0
    ggbld <- ggplot2::ggplot_build(gg)
    xticks = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.major_source # needed for <3.3.0
    xticks.minor = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.minor_source # needed for <3.3.0
    xlims = ggbld$layout$panel_params[[1]]$x.range

    # add positions of annotation labels
    xticks = c(xlims[1], xticks, xlims[2])
    # get breaks labels
    xtlabs = ggbld$layout$coord$labels(ggbld$layout$panel_params)[[1]]$x.labels # needed for <3.3.0
  }

  # align with \n
  txt <- c(left, right)
  txt = paste0("\n", txt)
  xtlabs = paste0(xtlabs, "\n")
  xtlabs = c(txt[1], xtlabs, txt[2])

  # return
  gg + ggplot2::scale_x_continuous(breaks = xticks, labels = xtlabs, minor_breaks = xticks.minor)

}


# extracts variables from a list of terms
mti_extract_variables <- function(lst) {
  # filter down only to the variables needed for plotting
  # need to parse x and ... list
  # browser()
  # q <- quos(...)
  vars = c()
  if (length(q) > 0) {
    vars <- lst %>% unlist() %>% lapply(function(x){x %>% all.vars()}) %>% unlist() %>% as.vector()
  }
  vars <- unique(vars)

  # return
  vars
}



#' Fixes the problem of exploding environments when saving ggplots to files
#'
#' Magic code by Mustafa (- JK)
#'
#' @param x ADD PARAM DESCRIPTION
#'
#' @return ADD RETURN DESCRIPTION
#'
#'
#' @author MB (JK)
#'
#' @import ggplot2
#'
#' @noRd
mti_fix_ggplot_env <- function(p) {
  # all the environment for quoted variables leads explosion of the object size
  # problem is also not that simple to just find those variables and clean up the
  # respective environment which I did, which did not solve the problem.
  # best solution so far is to wrap plot, get rid of everything else
  local(
    # transformartion of images into blank panel
    # this is updated
    ggplot2::ggplot(data.frame(x = 0:1, y = 0:1), ggplot2::aes_(x = ~x, y = ~y)) +
      ggplot2::geom_blank() +
      ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      ggplot2::annotation_custom(gg_grob, xmin = 0, xmax = 1 , ymin = 0, ymax = 1) +
      ggplot2::theme_void(),
    as.environment(list(gg_grob =ggplotGrob(p))) %>% {parent.env(.)=.GlobalEnv;.})
}


#' ADD TITLE
#'
#' ADD DESCRIPTION
#'
#' @param x ADD PARAM DESCRIPTION
#'
#' @return ADD RETURN DESCRIPTION
#'
#' @author whoWroteIt?
#' @noRd
fixorder = function(x){o= unique(as.character(x)); gdata::reorder.factor(x, new.order=o)} # fix order of a factor


#'
#'
#'
#'
#' @param Ds Concatenated dataframe returned by mti_format_se_samplewise
#' @param sample_filter term which samples to filter to first
#'
#' @return logical vector of samples to keep
#'
#' @author JK, JZ, KC
#' @noRd
mti_filter_samples <- function(Ds, filter_q, num_samp){

  Ds <- Ds %>%
    dplyr::mutate(tmpsamplenum = 1:nrow(Ds)) %>%
    dplyr::filter(!!filter_q) %>%
    droplevels()
  # message("filter metabolites: ", metab_filter_q, " [", nrow(stat), " remaining]")
  # did we leave 0 rows?
  if (nrow(Ds)==0) stop("Filtering left 0 rows")
  if (nrow(Ds)==num_samp) MetaboTools:::mti_logwarning('filtering did not filter out any samples')

  # store used samples
  samples.used <- rep(F, num_samp)
  samples.used[Ds$tmpsamplenum] <- T

  # return
  samples.used

}

#' Calculate evaluation measures for classifier predictions
#'
#' Calculates the following evaluation measures: sensitivity (sens), specificity (spec),
#'  accuracy (acc), positive-predictive value (ppv), and F1-measure (F1)
#'
#' @param tpfn a list of the four confusion matrix outcomes: True Positive (TP), False Positive (FP), False Negative (FN), and
#'  True Negative (TN)
#'
#' @return list of evaluation measures
#' @noRd
mti_measures <- function(tpfn) {
  # construct sensitivity, specificity, accuracy, PPV and F1 score; return as list
  x = list(
    sens = tpfn$TP/(tpfn$TP+tpfn$FN),
    spec = tpfn$TN/(tpfn$TN+tpfn$FP),
    acc = (tpfn$TP+tpfn$TN)/(tpfn$TP+tpfn$FN+tpfn$FP+tpfn$TN),
    PPV = tpfn$TP/(tpfn$TP+tpfn$FP)
  )
  x$F1 = 2*(x$sens*x$spec)/(x$sens+x$spec)
  x
}

#' Calculate Confusion Matrix Outcomes
#'
#' Calculates the following four confusion matrix outcomes: True Positive (TP), False Positive (FP), False Negative (FN), and
#'  False Positive (FP)
#'
#' @param trueclass a boolean vector of class labels
#' @param predclass a boolean vector of predicted classes
#'
#' @return list of confusion matrix outcomes
#' @noRd
mti_TPFN <- function(trueclass, predclass) {
  list(
    TP = sum(trueclass & predclass),
    FP = sum(!trueclass & predclass),
    FN = sum(trueclass & !predclass),
    TN = sum(!trueclass & !predclass)
  )
}

#' Get Evaluation Measures List
#'
#' Get a list of evaluation measures given a vector of predicted class probabilities and class labels. The following
#' evaluation measures are included in the list: sensitivity (sens), specificity (spec), accuracy (acc), positive-predictive
#' value (ppv), and F1-measure (F1).
#'
#' @param trueclass a boolean vector of class labels
#' @param predprob a vector of predicted class probabilities
#'
#' @return
#' @noRd
mti_get_measures_list <- function(trueclass, predprob) {
  # initialize arrays
  sensvals = c()
  specvals = c()
  ppvvals = c()
  accvals = c()
  f1vals = c()
  # loop over all predprob values as cutoff
  predprob.sorted = sort(predprob)
  # trick: add -Inf to make the curve work
  predprob.sorted = c(-Inf, predprob.sorted)
  # loop over all values
  for (i in 1:length(predprob.sorted)) { # could also be done via sapply()
    # do the cut
    predclass = predprob > predprob.sorted[i]
    # quality measures
    m = mti_measures(mti_TPFN(trueclass, predclass))
    sensvals[i] = m$sens
    specvals[i] = m$spec
    ppvvals[i] = m$PPV
    accvals[i] = m$acc
    f1vals[i] = m$F1
  }
  # calculate AUC under ROC curve
  # we have to take the negative, because we built the curve the backwards
  AUC = -pracma::trapz(1-specvals,sensvals)

  predprob.sorted[1] <- 0
  # return list
  list(sensvals=sensvals, specvals=specvals, AUC=AUC, ppvvals=ppvvals, thresholds=predprob.sorted, accvals=accvals, f1vals=f1vals)
}

#' Check if Data is Logged
#'
#' Check to see if either mt_pre_trans_log or mt_flag_logged was called in the pipeline. If mt_pre_trans_log called multiple
#' times or logging reversed, throw a warning and set is_logged to FALSE.
#'
#' @param D SummarizedExperiment
#'
#' @return is_logged (boolean)
#' @noRd
mti_check_is_logged <- function(D){

  is_logged <- FALSE

  # get names of the functions called
  called_functions <- names(metadata(D)$results)

  # if data is logged more than once, throw warning and return is_logged = FALSE
  num_times_logged <- called_functions %>% startsWith("pre_trans_log") %>% sum()
  if(num_times_logged > 1){
    warning("Data logged multiple times! Cannot determine log status. is_logged will be set to FALSE.")
    return(is_logged)
  }

  # if logging is reversed, throw warning and return is_logged = FALSE
  if(num_times_logged==1){
    log_index <- called_functions %>% startsWith("pre_trans_log") %>% which()
    exp_after_logged <- called_functions[(first_log_index+1):length(called_functions)] %>% startsWith("pre_trans_exp") %>% any()
    if(exp_after_logged){
      warning("Logging of data has been reversed! Cannot determine log stats. is_logged will be set to FALSE.")
      return(is_logged)
    }else{
      is_logged <- TRUE
    }
  }

  if(any(startsWith(called_functions, "flag_logged"))){
    # if pre_trans_log called AND flag_logged called, return is_logged as False
    if(is_logged){
      is_logged <- FALSE
      return(is_logged)
    }else{
      is_logged <- TRUE
    }
  }

  is_logged

}
