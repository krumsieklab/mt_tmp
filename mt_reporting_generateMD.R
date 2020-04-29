
#' Markdown-based report generator
#' 
#' Generates a fully automated report version of an entire (linear) pipeline.  
#' 
#' Will generate RMD markdown document that can be adapted and then knitted to HTML.
#'
#' 
#' ### Markdown generation
#' EXAMPLE COMING SOON
#' 
#' 
#' @author JK
mt_reporting_generateMD <- function(
  D,                     # SE
  outfile = 'MT.RMD',    # output file to write to
  readfrom = 'mt.rds',   # R data file data will be loaded and is supposed to contain SummarizedExperiment "D"
  title = 'RMD output',  # title of document
  firstheading='Output', # name of first heading
  use.plotly=F,          # output interactive plotly plots?
  output.calls=F,        # output full information about all parameters of each function call?
  start.after=NA         # UUID of pipeline step AFTER which to start (default: none, i.e. output entire pipeline) 
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
  h <- file(outfile, open='wt')
  
  #### markdown header and first heading
  out(glue('
---
title: {title}
output: 
  html_document:
    toc: true
    toc_float: TRUE
---

    '))
  
  #### determine where to start from (first pipeline step, or after one)
  start.from = 1 # by default
  if (!is.na(start.after)) {
    # extract all UUIDs
    allids <- D %>% metadata() %>% .$results %>% map("uuid") %>% unlist()
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
    writechunk('# load libraries\nmt.quickload()\nlibrary("DT")')  
  } else  {
    writechunk('# load libraries\nmt.quickload()\nlibrary("DT")\nlibrary("plotly")')  
  }
  #### chunk that loads data
  writechunk(glue('# load data\nload("{readfrom}")\nr <- metadata(D)$results'))
  
  #### start of output
  out(glue('# {firstheading}'))
  
  #### loop over results
  lvl=2 # current level of results = 2
  
  # loop over results
  r <- metadata(D)$results
  
  for (i in start.from:length(r)) {
    # to ignore?
    if (r[[i]]$fun[2]!="void") { # ignore void
      
      # reporting step?
      if (r[[i]]$fun[1]!="reporting") {
        # not reporting, actual pipeline step
        
        ## header
        out(glue('{strrep("#",lvl)} {r[[i]]$fun %>% paste(collapse="_")}'))
        
        ## detailed arguments?
        if (output.calls) {
          L <- r[[i]]$args 
          out("*Function arguments:*<br/>")
          out(names(L) %>% lapply(function(x){sprintf("%s=%s",x, toString(L[[x]]))}) %>% paste0(collapse = "<br/>"))
          out("")
        }
        
        ## log text
        out(glue('*Log text:*<br/>{r[[i]]$logtxt}\n\n'))
        
        ## plot?
        if (r[[i]]$fun[1]=="plots") {
          # plot
          if (!use.plotly) {
            writechunk( glue("r[[{i}]]$output"))
          } else {
            writechunk( glue("
plotlist = r[[{i}]]$output %>% lapply(ggplotly)
htmltools::tagList(setNames(plotlist, NULL))
                             "), params='results="show"')
          }
          # }
          
        }
        
        ## statistical result table?
        if (r[[i]]$fun[1]=="stats") {
          # write out datatable
          writechunk(glue('
# extract result table
df<-r[[{i}]]$output$table
# add metabolite names
rd <- rowData(D)
df <- cbind(name=as.data.frame(rd)$name[match(df$var, rownames(rd))], df)
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
          out(glue('{strrep("#",r[[i]]$output$lvl)} {r[[i]]$output$title}'))
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