#' Generates markdown-based HTML output for Non-Linear Pipelines from a list of SummarizedExperiment Objects
#'
#' Shortcut for mt_reporting_generateMD_nonLinear with subsequent knitting and clean-up.
#'
#' @param pipelines A list of SummarizedExperiment objects
#' @param outfile # output HTML file name
#' @param title Title of RMD document
#' @param output.calls Output detailed info on function calls? default: F (passed through to mt_reporting_generateMD)
#' @param keep.tmp Keep temporary files? (default: No). Can be used to manually edit RMD afterwards.
#'
#' @author JK, KC
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_reporting_html_nonLinear <- function(
  pipelines,
  outfile,
  title = 'RMD output',
  output.calls=F,
  keep.tmp=F
) {
  
  # validate argument
  ## pipelines
  stopifnot("list" %in% class(pipelines))
  check_SE <- sapply(pipelines, function(p){"SummarizedExperiment" %in% class(p)})
  if(any(check_SE==FALSE)){
    stop("Argument pipelines must be a list of SummarizedExperiment objects.")
  }
  ## outfile
  if(missing(outfile)){
    stop("Argument outfile must be provided.")
  }
  
  
  res <- MTResultCollector$new()
  res$addMultiple(pipelines)
  
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
    output_calls = output.calls)
  # save temp file that will be input for the RMD
  save(D, file=rdsfile)
  # knit
  rmarkdown::render(rmdfile)
  # rename to correct name
  file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), outfile)
  # clean up
  if (!keep.tmp) {
    file.remove(rmdfile)
    file.remove(rdsfile)
  }
  
}
