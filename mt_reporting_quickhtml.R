#' Generates markdown-based HTML output from SummarizedExperiment
#' 
#' Shortcut for mt_reporting_generateMD with subsequent knitting and clean-up.
#'
#' @param D Summarized experiment input
#' @param outfile # output HTML file name
#'
#' @author JK
#' 
mt_reporting_quickhtml <- function(D, outfile, title = 'RMD output') {
  
  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))  
  
  # define file names
  rmdfile <- "tmpmd.RMD"
  rdsfile <- "tmpmd.rds"
  # generate RMD
  D %>% mt_reporting_generateMD(outfile = rmdfile, readfrom = rdsfile, title = title)
  # save temp file that will be input for the RMD
  save(D, file=rdsfile)
  # knit
  rmarkdown::render(rmdfile)
  # rename to correct name
  file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), outfile)
  # clean up
  file.remove(rmdfile)
  file.remove(rdsfile)
  
}
