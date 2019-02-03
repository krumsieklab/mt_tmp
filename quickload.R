# Temporary solution until MT becomes a package.
# Source this file to load all functions

library(tidyverse)

if(!exists("codes.makepath"))
  codes.makepath <- function(a)paste0("./", a)

l <- list.files(path=codes.makepath("packages/metabotools"),pattern="*.R$",full.names=T)
l <- l[!grepl('*examples*|*quickload*',l)]
suppressPackageStartupMessages({
  walk(l, function(x){source(x,echo=F,verbose=F)})
})
