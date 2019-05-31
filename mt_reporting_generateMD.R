
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
  firstheading='Output'  # name of first heading
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
  
  
  
  #### global chunk options
  writechunk("# default chunk options\nknitr::opts_chunk$set(warning=F,echo=F,results='hide',message=F)", params = "echo=F")  
  
  #### chunk that loads libraries
  writechunk('# load libraries\nmt.quickload()')  
  #### chunk that loads data
  writechunk(glue('# load data\nload("{readfrom}")\nr <- metadata(D)$results'))
  
  #### start of output
  out(glue('# {firstheading}'))
  
  #### loop over results
  lvl=2 # current level of results = 2
  
  # loop over results
  r <- metadata(D)$results
  
  for (i in 1:length(r)) {
    # reporting step?
    if (r[[i]]$fun[1]!="reporting") {
      # not reporting, actual pipeline step
      
      ## header
      out(glue('{strrep("#",lvl)} {r[[i]]$fun %>% paste(collapse="_")}'))
      
      ## log text
      out(glue('{r[[i]]$logtxt}\n'))
      
      ## plot?
      if (r[[i]]$fun[1]=="plots") {
        # # test code for subchunkify
        # if (all(all.equal(r[[i]]$fun, c("plots","qc","dilutionplot")) == T)) {
        #   # quot norm plot, try different size
        #   writechunk( glue("r[[{i}]]$output %>% lapply(subchunkify, fig_height=30, fig_width=30 )"))
        # } else {
        # any other plot
        writechunk( glue("r[[{i}]]$output %>% lapply(plot)"))
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
df <- cbind(name=rd$name[match(df$var, rownames(rd))], df)
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
    }
  }
  
  # clean up
  close(h)
  
}