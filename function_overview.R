# Generates an overview of all functions including they Roxygen titles (first line).
# Also determines which example scripts use which function.
# All done via regular expression matching.

#### initialize ----
library(tidyverse)
library(stringr)
library(openxlsx)
setwd(codes.makepath("MT"))

# read file into string
slurp <- function(file){readChar(file, file.info(file)$size)}

#### extract information

## Parse out all functions and their roxygen titles
titles <-
  # list of all code files
  list.files(path = "MetaboTools/R", pattern = "mt.*.R", full.names = T) %>%
  # no internal files
  {x<-.; x<-x[!grepl("internal",x)]; x} %>%
  # extract all Roxygen titles (first line, everything after #')
  lapply(function(file){
    # for each function in the file, extract the whole block from first Roxygen line #' to function name
    chunks <- slurp(file) %>% str_extract_all(regex("^#'.*?mt_.*?<-.*?function", multiline=TRUE, dotall = T)) %>% .[[1]]
    # for each of those chunks, get the first Roxygen comment (= the title) and the function name
    chunks %>% lapply(function(chunk){
      data.frame(
        fun = chunk %>% str_match("(mt_.*?)<-") %>% .[1,2] %>% trimws(),
        title = chunk %>% str_match("#'(.*)") %>% .[1,2] %>% trimws()
      )
    }) %>% do.call(rbind,.)
  })  %>% do.call(rbind,.)
# get rid of internal functions
# remove path and ".R" in all names
names(titles) %<>% basename() %>% gsub("\\.R", "", .)

## find which examples call which functions
excalls  <-
  # list all example code files
  list.files(path = "examples", pattern = "*.R", full.names = T) %>%
  # extract all functios called by this example function
  sapply(function(file){slurp(file) %>% str_match_all("(mt_.*?)\\(") %>% .[[1]] %>% .[,2]})
# remove path and ".R" in all names
names(excalls) %<>% basename() %>% gsub("\\.R", "", .)

## build data frame that contains function names, titles, and which example files these the functions are called in
names(titles) %>% lapply(function(fun){
  cbind(
    # function name and title
    data.frame(fun=fun, title=titles[[fun]]),
    # matching function name in each example list
    excalls %>% sapply(function(x){ifelse(!is.na(any(match(x,fun))),"X","")}) %>% as.list()
  )
}) %>% do.call(rbind, .) %>%
  # write out to Excel
  openxlsx::write.xlsx(file="function_overview.xlsx")

