#' Generates markdown-based HTML output from SummarizedExperiment
#' 
#' Shortcut for mt_reporting_generateMD with subsequent knitting and clean-up.
#'
#' @param D Summarized experiment input
#' @param outfile # output HTML file name
#' @param output.calls Output detailed info on function calls? default: F (passed through to mt_reporting_generateMD)
#' @param start.after UUID of pipeline step AFTER which to start (default: none, i.e. output entire pipeline) (passed through to mt_reporting_generateMD)
#'
#' @author JK
#' 
mt_reporting_quickhtml <- function(
  D, 
  outfile, 
  title = 'RMD output',
  output.calls=F,
  start.after=NA
  ) {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))  
  
  # define file names
  rmdfile <- "tmpmd.RMD"
  rdsfile <- "tmpmd.rds"
  # generate RMD
  D %>% mt_reporting_generateMD(outfile = rmdfile, readfrom = rdsfile, title = title, output.calls = output.calls, start.after=start.after)
  # save temp file that will be input for the RMD
  save(D, file=rdsfile)
  # knit
  rmarkdown::render(rmdfile)
  # rename to correct name
  file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), outfile)
  # clean up
  file.remove(rmdfile)
  file.remove(rdsfile)
  
  # return document, in case pipeline is supposed to keep running
  D
}
