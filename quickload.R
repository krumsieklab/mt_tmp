# Temporary solution until MT becomes a package.
# Source this file to load all functions

library(tidyverse)

l <- list.files(path=codes.makepath("MT"),pattern="*.R$",full.names=T)
l <- l[!grepl('*examples*|*quickload*|*ignore*',l)]
suppressPackageStartupMessages({
  walk(l, function(x){source(x,echo=F,verbose=F)})
})

print("MT sourced.")
