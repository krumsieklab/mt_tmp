# # MetaboTools
#
# Helper functions.
#
# last update: 2019-3-15
# authors: JK,MB
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
  d <- D %>% mti_format_se_samplewise()
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




