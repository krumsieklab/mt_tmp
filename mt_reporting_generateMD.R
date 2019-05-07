#### Markdown-based report generator



mt_reporting_generateMD <- function(D, firstheading='Output', directknit=T) {
  
  #### define inside functions
  
  # write sprintf style to the console
  writecons <- function(txt)cat(paste0(txt,"\n"))
  # write sprintf style to the currently open file
  writefile <- function(txt)writeLines(txt, stdout())
  # write either to console or file, depending on directknit 
  writeswitch <- function(txt)if(directknit){writecons(txt)}else{writefile(txt)}
  # generate an R chunk
  writechunk <- function(code) {
    writefile("```{r}")
    writefile(code)
    writefile("```")
  }
  # either execute command or write it out as chunk
  execswitch <- function(cmd)if(directknit){eval(parse(text=cmd))}else{writechunk(cmd)}

  
  #### start actual output
  
  # first heading
  writeswitch(glue('# {firstheading}'))
  lvl=2 # current level of results = 2
  
  # loop over results
  r <- metadata(D)$results
  
  for (i in 1:length(r)) {
    # reporting step?
    if (r[[i]]$fun[1]!="reporting") {
      # not reporting, actual pipeline step
      
      ## header
      writeswitch(glue('{strrep("#",lvl)} {r[[i]]$fun %>% paste(collapse="_")}'))
      
      ## log text
      writeswitch(glue('{r[[i]]$logtxt}\n'))
      
      ## plot?
      if (r[[i]]$fun[1]=="plots") {
        execswitch( glue("r[[{i}]]$output %>% lapply(plot)"))
      }
      
      ## statistical result table?
      if (r[[i]]$fun[1]=="stats") {
        execswitch(glue(' r[[{i}]]$output %>% .$table %>% knitr::kable(format="html") %>% cat()'))
      }
      
      # footer
      writeswitch("\n")
      
    } else {
      
      # special reporting step
      if (r[[i]]$fun[2]=="heading") {
        # add extra heading
        writeswitch(glue('{strrep("#",r[[i]]$output$lvl)} {r[[i]]$output$title}'))
        # result level is this heading +1
        lvl = r[[i]]$output$lvl + 1
      }
    }
  }
  
}

# D %>% mt_reporting_generateMD(directknit = F)

# knitr::kable(output, "html", booktabs = TRUE, longtable = TRUE, caption = "Test") %>%
# kableExtra::kable_styling(latex_options = c("hold_position", "repeat_header"))
# %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
