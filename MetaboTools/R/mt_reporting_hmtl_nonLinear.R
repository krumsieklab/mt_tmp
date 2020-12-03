# QUESTIONS: DELETE ONCE ANSWERED
# 1. Should this return something? mt_reporting_html returns D, but this works on a LIST of objects, so can't be part of a pipeline

# FIXES
# 1. Need to update description when integrate mt_reporting_generateMD
# 2. Add @examples

#' Generates markdown-based HTML output for non-linear pipelines from a list of SummarizedExperiment objects
#'
#' Shortcut for mt_reporting_generateMD_nonLinear with subsequent knitting and clean-up.
#'
#' @param D_list A list of \code{SummarizedExperiment} objects.
#' @param outfile Output HTML filename.
#' @param title Title of RMD document. Default: 'Non-Linear RMD output'.
#' @param output_calls Output detailed info on function calls? Default: F.
#' @param keep_tmp Keep temporary files? Can be used to manually edit RMD afterwards. Default: F.
#'
#' @return Nothing. This does not pass through \code{SummarizedExperiment} objects.
#'
#' @author JK, KC
#'
#' @export
mt_reporting_html_nonlinear <- function(D_list,
                                        outfile,
                                        title = 'Non-Linear RMD output',
                                        output_calls=F,
                                        keep_tmp=F
) {
  
  # validate argument
  ## D_list
  stopifnot("list" %in% class(D_list))
  check_SE <- sapply(D_list, function(p){"SummarizedExperiment" %in% class(p)})
  if(any(check_SE==FALSE)){
    stop("Argument D_list must be a list of SummarizedExperiment objects.")
  }
  ## outfile
  if(missing(outfile)){
    stop("Argument outfile must be provided.")
  }
  
  
  res <- MTResultCollector$new()
  res$addMultiple(D_list)

  # throw error if multiple roots (i.e. disjointed D_list)
  if(length(res$graph_roots()) > 1){
    stop("D_list can not be disjointed! A common root must exist for all D_list!")
  }

  # unique string
  ustr <- uuid::UUIDgenerate()
  
  # define file names
  rmdfile <- sprintf("tmp_%s.RMD", ustr)
  rdsfile <-  sprintf("tmp_%s.rds", ustr)
  # generate RMD
  res %>% mt_reporting_generateMD_nonLinear(
    outfile = rmdfile,
    read_from = rdsfile,
    title = title,
    output_calls = output_calls)
  # save temp file that will be input for the RMD
  save(D, file=rdsfile)
  # knit
  rmarkdown::render(rmdfile)
  # rename to correct name
  file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), outfile)
  # clean up
  if (!keep_tmp) {
    file.remove(rmdfile)
    file.remove(rdsfile)
  }
  
}
