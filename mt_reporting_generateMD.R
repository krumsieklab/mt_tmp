#### Markdown-based report generator

cats <- function(...)cat(paste0(sprintf(...),"\n"))

mt_reporting_generateMD <- function(D, firstheading='Output') {
  
  
  # first heading
  cats("# %s", firstheading)
  lvl=2 # current level of results = 2
  
  # loop over results
  r <- metadata(D)$results
  
  for (i in 1:length(r)) {
    # reporting step?
    if (r[[i]]$fun[1]!="reporting") {
      # not reporting, actual pipeline step
      # header
      cats("%s %s", strrep("#",lvl), r[[i]]$fun %>% paste(collapse="_"))
      # log text
      cats("%s\n", r[[i]]$logtxt)
      # plot?
      if (r[[i]]$fun[1]=="plots") {
        r[[i]]$output %>% lapply(plot)
      }
      # statistical result table?
      if (r[[i]]$fun[1]=="stats") {
        cats(
          r[[i]]$output %>% 
            .$table %>% 
            knitr::kable()
        ) 
      }
      
      # footer
      cats("\n")
      
    } else {
      # special reporting step
      if (r[[i]]$fun[2]=="heading") {
        # add extra heading
        cats("%s %s", strrep("#",r[[i]]$output$lvl), r[[i]]$output$title)
        # result level is this heading +1
        lvl = r[[i]]$output$lvl + 1
      }
    }
  }
  
}


# knitr::kable(output, "html", booktabs = TRUE, longtable = TRUE, caption = "Test") %>%
# kableExtra::kable_styling(latex_options = c("hold_position", "repeat_header"))
# %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
