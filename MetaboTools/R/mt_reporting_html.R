# QUESTIONS: DELETECE ONCE ANSWERED
# 1. Can (experimental) be deleted?

# FIXES
# 1. Need to update description when integrate mt_reporting_generateMD
# 2. Add @examples

#' Generates markdown-based HTML output from SummarizedExperiment
#'
#' Shortcut for mt_reporting_generateMD with subsequent knitting and clean-up.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Output HTML filename.
#' @param title Title of RMD document. Default: 'RMD output'.
#' @param output_calls Output detailed info on function calls? Default: F.
#' @param number_sections Number sections and sub-sections? Default: F.
#' @param start_after UUID of pipeline step AFTER which to start. Default: NA, i.e. output entire pipeline. (passed through to mt_reporting_generateMD)
#' @param use_plotly Output interactive plotly plots? Default: F. (experimental)
#' @param keep_tmp Keep temporary files? Can be used to manually edit RMD afterwards. Default: F.
#'
#' @return Does not change the \code{SummarizedExperiment} object.
#'
#' @author JK
#'
#' @export
mt_reporting_html <- function(D,
                              file,
                              title = 'RMD output',
                              output_calls=F,
                              number_sections=F,
                              start_after=NA,
                              use_plotly=F,  # EXPERIMENTAL
                              keep_tmp=F
) {

  # validate argument
  stopifnot("SummarizedExperiment" %in% class(D))

  # unique string
  ustr <- uuid::UUIDgenerate()

  # define file names
  rmdfile <- sprintf("tmp_%s.RMD", ustr)
  rdsfile <-  sprintf("tmp_%s.rds", ustr) # only used of keep_tmp==T

  # generate RMD
  D %>% mt_reporting_generateMD(
    file = rmdfile,
    readfrom = rdsfile,
    title = title,
    output.calls = output_calls,
    number.sections = number_sections,
    start.after=start_after,
    use.plotly = use_plotly)

  # save temp file that will be input for the RMD?
  if (keep_tmp) save(D, file=rdsfile)
  # knit
  rmarkdown::render(rmdfile, params=list(D=D))
  # rename to correct name
  file.rename(paste0(tools::file_path_sans_ext(rmdfile),'.html'), file)

  # clean up
  if (!keep_tmp) {
    # no temp files left behind
    file.remove(rmdfile)
    # .rds does not need to be deleted because it was never generated
  } else {
    # keep temp files, move them to their final destination name
    file.rename(rmdfile, paste0(tools::file_path_sans_ext(file),'.rmd'))
    file.rename(rdsfile, paste0(tools::file_path_sans_ext(file),'.rds'))
  }

  # return document, in case pipeline is supposed to keep running
  D
}
