#' Markdown-based report generator
#'
#' Generates a fully automated report version of an entire (linear) pipeline.
#'
#' Will generate RMD markdown document that can be adapted and then knitted to HTML.
#'
#' @param D \code{SummarizedExperiment} input
#' @param file File to be generated
#' @param readfrom Name of R data file data will be loaded and is supposed to contain SummarizedExperiment "D". Will not actually be loaded in this function, but while knitting the RMD later.
#' @param title Title of RMD document
#' @param firstheading Name of first heading
#' @param use.plotly Output interactive plotly plots? (experimental)
#' @param output.calls Output full information about all parameters of each function call into RMD?
#' @param number.sections Number sections and sub-sections? (default: F)
#' @param start.after UUID of pipeline step AFTER which to start (default: none, i.e. output entire pipeline)
#'
#' @returns nothing
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#'
#' @author JK
#'
#' @noRd

mt_reporting_generateMD <- function(
  D,
  file = 'MT.RMD',
  readfrom = 'mt.rds',
  title = 'RMD output',
  firstheading='Output',
  use.plotly=F,
  output.calls=F,
  number.sections=F,
  start.after=NA
) {


  #### helper functions
  out <- function(str){writeLines(str, h)}

  writechunk <- function(code, params='') {
    if (nchar(params)>0) params=paste0(' ', params)
    out(sprintf("```{r%s}", params))
    out(code)
    out("```\n")
  }


  #### initialize output file
  h <- file(file, open='wt')

  #### markdown header and first heading
  out(glue::glue('
---
title: {title}
output:
  html_document:
    toc: true
    toc_float: TRUE
    {if(number.sections){"number_sections: true"}else{""}}
params:
  D: NA
---

    '))

  #### determine where to start from (first pipeline step, or after one)
  start.from = 1 # by default
  if (!is.na(start.after)) {
    # extract all UUIDs
    allids <- D %>% metadata() %>% .$results %>% purrr::map("uuid") %>% unlist()

    # if tag name provided, get uuid
    tag_name_list <- MetaboTools:::mtm_res_get_path(D, c("reporting", "tag")) %>% purrr::map("output") %>% unlist()
    tag_name_idx <- match(start.after, tag_name_list)
    if(!is.na(tag_name_idx)){
      start.after <- tag_name_list[tag_name_idx] %>% names() %>% gsub("reporting_tag.", "", .)
    }

    # find the one to start after
    start.from <- which(allids==start.after)
    # error-check
    if (length(start.from)==0)
      stop(sprintf("Could not find pipeline step to start after: '%s'", start.after))
    if (D %>% metadata() %>% .$results %>% length() == start.from)
      stop(sprintf("Cannot start after pipeline step '%s', because it's the last entry of the pipeline", start.after))
    start.from <- start.from + 1
  }



  #### global chunk options
  writechunk("# default chunk options\nknitr::opts_chunk$set(warning=F,echo=F,results='hide',message=F)", params = "echo=F")

  #### chunk that loads libraries
  if (!use.plotly) {
    writechunk('# load libraries\nlibrary(MetaboTools)\n')
  } else  {
    writechunk('# load libraries\nlibrary(MetaboTools)\nlibrary("plotly")')
  }
  #### chunk that assigns list r from data
  writechunk(glue::glue('r <- metadata(D)$results'))

  #### start of output
  out(glue::glue('# {firstheading}'))

  #### loop over results
  lvl=2 # current level of results = 2

  # loop over results
  r <- metadata(D)$results

  for (i in start.from:length(r)) {
    # to ignore?
    if (length(r[[i]]$fun)<2 || r[[i]]$fun[2]!="void") { # ignore void

      # reporting step?
      if (r[[i]]$fun[1]!="reporting") {
        # not reporting, actual pipeline step

        ## header
        out(glue::glue('{strrep("#",lvl)} {r[[i]]$fun %>% paste(collapse="_")}'))

        ## detailed arguments?
        if (output.calls) {
          L <- r[[i]]$args
          out("*Function arguments:*<br/>")
          out(names(L) %>% lapply(function(x){sprintf("%s=%s",x, toString(L[[x]]))}) %>% paste0(collapse = "<br/>"))
          out("")
        }

        ## log text
        out(glue::glue('*Log text:*<br/>{r[[i]]$logtxt}\n\n'))

        ## plot?
        if (r[[i]]$fun[1]=="plots") {

          # special parameters?
          extraparams <- ""
          if (r[[i]]$fun[2]=="statsbarplot") {
            # dynamic height
            if(r[[i]]$output2$nr!=0){
              # set plot height
              height <- (62+(r[[i]]$output2$nr*r[[i]]$output2$npanrow*23.7143)+84)/2304*12 # manually curated using pixel measurements on example
              width <- 5+3*r[[i]]$output2$npancol
            } else{
              # if empty plot, set height to 3
              height <- 3
              width <- 8
            }
            extraparams <- sprintf(",fig.width=%f,fig.height=%f", width, height)
          } else if(length(r[[i]]$fun)>=3){
            if(r[[i]]$fun[2]=="boxplot"&r[[i]]$fun[3]=="scatter"){
              # dynamic height
              if(!is.null(r[[i]]$output2)){
                # set plot height
                height <- (23+(r[[i]]$output2*160.65)+33)/2304*32 # manually curated using pixel measurements on example
                width <- 7
              } else{
                # if empty plot, set height to 3
                height <- 3
                width <- 7
              }
              extraparams <- sprintf(",fig.width=%f,fig.height=%f", width, height)
            }
          }



          # plot
          if (!use.plotly) {
            writechunk( glue::glue("r[[{i}]]$output"), params = extraparams)
          } else {
            writechunk( glue::glue("
plotlist = r[[{i}]]$output %>% lapply(ggplotly)
htmltools::tagList(setNames(plotlist, NULL))
                             "), params=paste0('results="show"', extraparams))
          }
          # }

        }

        ## statistical result table?
        if (r[[i]]$fun[1]=="stats") {
          # write warning if df too large
          if(nrow(r[[i]]$output$table) > 1000){
            out(glue::glue('WARNING: Large data frame ({nrow(r[[i]]$output$table)} rows). Displaying first
                           1000 rows.'))
          }
          # write out datatable
          writechunk(glue::glue('
# extract result table
df<-r[[{i}]]$output$table
# add metabolite names
rd <- rowData(D)
df <- cbind(name=as.data.frame(rd)$name[match(df$var, rownames(rd))], df) %>%
  dplyr::arrange(p.value)
# subset large data frames
if(nrow(df) > 1000) df <- df[1:1000, ]
# output
DT::datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 20, lengthMenu = c(10*(2^(0:3)), nrow(df)), autoWidth = TRUE, width = 1200, dom = "Bitlrp", buttons = c("copy", "csv", "excel", "pdf", "print")), class = "cell-border stripe", extensions = "Buttons")  %>% DT::formatStyle(columns = c(1:ncol(df)), fontSize = "80%", target= "row", lineHeight="80%")'),
                     params = "results='asis'")

        }

        # empty line as spacer
        out("")


      } else {

        # special reporting step
        if (r[[i]]$fun[2]=="heading") {
          # add extra heading
          out(glue::glue('{strrep("#",r[[i]]$output$lvl)} {r[[i]]$output$title}'))
          # result level is this heading +1
          lvl = r[[i]]$output$lvl + 1
        }

        # special reporting step
        if (r[[i]]$fun[2]=="text") {
          # enter text
          out(r[[i]]$output$text)
          out("")
        }
      }
    }
  }

  # clean up
  close(h)

}
