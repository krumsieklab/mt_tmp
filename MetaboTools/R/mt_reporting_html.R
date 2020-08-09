#' Generates markdown-based HTML output from SummarizedExperiment
#'
#' Shortcut for mt_reporting_generateMD with subsequent knitting and clean-up.
#'
#' @param D Summarized experiment input
#' @param outfile # output HTML file name
#' @param title Title of RMD document
#' @param output.calls Output detailed info on function calls? default: F (passed through to mt_reporting_generateMD)
#' @param number.sections Number sections and sub-sections? (default: F)
#' @param start.after UUID of pipeline step AFTER which to start (default: none, i.e. output entire pipeline) (passed through to mt_reporting_generateMD)
#' @param use.plotly Output interactive plotly plots? (experimental)
#' @param keep.tmp Keep temporary files? (default: No). Can be used to manually edit RMD afterwards.
#'
#' @author JK
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @export

mt_reporting_html <- function(
  D,
  outfile,
  title = 'RMD output',
  output.calls=F,
  number.sections=F,
  start.after=NA,
  use.plotly=F,  # EXPERIMENTAL
  keep.tmp=F
) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # unique string
  ustr <- uuid::UUIDgenerate()

  # define file names
  rmdfile <- sprintf("tmp_%s.RMD", ustr)
  if(keep.tmp) rdsfile <-  sprintf("tmp_%s.rds", ustr)
  # generate RMD
  D %>% mt_reporting_generateMD(
    outfile = rmdfile,
    readfrom = rdsfile,
    title = title,
    output.calls = output.calls,
    number.sections = number.sections,
    start.after=start.after,
    use.plotly = use.plotly)
  # save temp file that will be input for the RMD
  if (keep.tmp) save(D, file=rdsfile)
  # knit
  rmarkdown::render(rmdfile, params=list(D=D))
  # rename to correct name
  file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), outfile)
  # clean up
  if (!keep.tmp) {
    file.remove(rmdfile)
  }

  # return document, in case pipeline is supposed to keep running
  D
}
