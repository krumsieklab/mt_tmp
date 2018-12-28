# MetaboTools
#
# Add pathway information
#
# Adds custom pathways to the already existing SummarizedExperiment
# data structure.
#
# last update: 2018-12-11
# authors: Parviz Gomari
#

# todo: document output

# dependencies
library(tidyverse)
library(magrittr)
library(graphite)
library(rpubchem)


# main
mt_add_pathways <- function(
  D,             # SummarizedExperiment input
  in_col,        # column name to use for pathway fetching
  out_col,       # name of the column to add pathway info
  pw_species,    # name of the species the data was measured in
  pw_name        # the name of the pathway database to use
) {
  
  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))
  
  
  ######################################################################################
  ## This will be relocated into another script where all ID combinations
  ## have been remapped to graphite databases.
  ######################################################################################
  # to see all available databasaes and species, use:
  # pathwayDatabases() # %>% filter(species == "hsapiens")
  pwdb <- pathways(species = pw_species, database = pw_name)
  
  # use given number of cores, this should be part of input
  options(Ncpus = 5)
  
  # convert from ChEBI IDs (given by graphite) to our lieblings ID
  con_id <- if_else(in_col == "KEGG", "KEGGCOMP", in_col) 
  pwdb %<>% convertIdentifiers(con_id)
  
  # graphite comes with as a list of pathways given a database
  # given this list, subselect the metabolite entries per pathway
  pwdb <- 
    mcmapply(function(pwname) {
      
      pw <- pwdb[[pwname]]
      pw %>% 
        edges(which = "metabolites") %>% 
        mutate(name = pwname,
               ID = pathwayId(pw))
      
    },
    names(pwdb),
    SIMPLIFY = F,
    mc.cleanup = T,
    mc.cores = 5)
  
  # convert the pathway list into a dataframe 
  pwdb %<>% 
    bind_rows() %>% 
    dplyr::select(src, dest, name, ID) %>% 
    gather(key = tmp, value = mappingID, -c(name, ID)) %>% 
    dplyr::select(mappingID, name, ID) %>% 
    distinct()
  
  # create a dataframe with only the ID and name
  pwdb_map <- select(pwdb, ID, name)
  
  # nest all the pathway IDs given our lieblings input IDs
  pwdb_reduced <- 
    pwdb %>% 
    group_by(mappingID) %>% 
    nest(ID, .key = IDs) %>% 
    mutate(IDs = as.list(unlist(IDs, recursive = FALSE)))
  
  ######################################################################################
  ## End of move
  ######################################################################################
  
  # have to do this to take a variable input name for the columns
  # used for joining in the next step
  left_index <- enquo(in_col)
  right_index <- "mappingID"
  by <-  set_names(quo_name(right_index), quo_name(left_index))
  
  # match the nested pathways to our lieblings IDs
  pw_col <- D %>% 
    rowData %>% 
    .[in_col] %>% 
    as.data.frame() %>% 
    left_join(pwdb_reduced, by = by) %>% 
    .$IDs
  
  # add the pathway IDs into D
  rowData(D)[[out_col]] <- pw_col
  
  # add pathway map to the metadata of D
  metadata(D)$pathways[[out_col]] <- pwdb_map
  
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('added pathway annotations using the %s database', pw_name)
    )
  
  D
}




# Example -----------------------------------------------------------------

D_alone <- 
  mt_files_load_metabolon(codes.makepath("packages/metabotools/sampledata.xlsx"), "OrigScale") %>% 
  mt_add_pathways(in_col = "KEGG", out_col = "humancyc_db", pw_species = "hsapiens", pw_name = "humancyc")



