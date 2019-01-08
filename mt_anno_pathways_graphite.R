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
  pw_name,       # the name of the pathway database to use
  n_cpus = 2     # number of cores to use for parallelizaion (used in convertIdentifiers)
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
  options(Ncpus = n_cpus)
  
  # convert from ChEBI IDs (given by graphite) to our lieblings ID
  con_id <- if_else(in_col == "KEGG", "KEGGCOMP", in_col) 
  pwdb %<>% convertIdentifiers(con_id)
  
  # graphite comes with as a list of pathways given a database
  # given this list, subselect the metabolite entries per pathway
  pwdb <- 
    mcmapply(function(pwname) {
      
      pw <- pwdb[[pwname]]
      pw %>% 
        graphite::edges(which = "metabolites") %>% 
        mutate(name = pwname,
               ID = pathwayId(pw))
      
    },
    names(pwdb),
    SIMPLIFY = F,
    mc.cleanup = T,
    mc.cores = n_cpus)
  
  # convert the pathway list into a dataframe 
  pwdb %<>% 
    bind_rows() %>% 
    dplyr::select(src, dest, name, ID) %>% 
    gather(key = tmp, value = mappingID, -c(name, ID)) %>% 
    dplyr::select(mappingID, name, ID) %>% 
    drop_na() %>% # ensure that there are no NAs in any rows
    distinct()
  
  # create a dataframe that enables the mapping between pathway
  # names and IDs. Included also are num_total, num_measured,
  # num_pw_total, and num_pw_measured (see below for further details)
  pwdb_summary <-
    pwdb %>% 
    mutate(
      # num_total - overall number of metabolites in entire database (this will 
      # be a redundant, repeating number, identical in every row… but it’s the 
      # easiest way to store it right now)
      num_total = n_distinct(mappingID),
      # num_measured - overall number of measured metabolites found in the DB
      num_measured = 
        mappingID %>% 
        intersect(rowData(D)[[in_col]]) %>% 
        length()
    ) %>% 
    group_by(ID) %>% 
    mutate(
      # num_pw_total - the total number of metabolites in that pathway 
      # (overall DB background)
      num_pw_total = n(),
      # num_pw_measured - the number of measured metabolites in that pathway 
      # (for the M type of analysis on the actual measured background).
      num_pw_measured = 
        mappingID %>% 
        intersect(rowData(D)[[in_col]]) %>% 
        length()
    ) %>% 
    ungroup() %>% 
    dplyr::select(-mappingID) %>% 
    distinct()
  
  
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
  
  # add pathway database to the metadata of D
  metadata(D)$pathways[[out_col]]$database <- pwdb
  
  # add pathway map to the metadata of D
  metadata(D)$pathways[[out_col]]$id_map <- pwdb_summary
  
  
  
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('added pathway annotations using the %s database', pw_name)
    )
  
  D
}



if (F) {
  # Example -----------------------------------------------------------------
  D_alone <- 
    mt_files_load_metabolon(codes.makepath("packages/metabotools/sampledata.xlsx"), "OrigScale") %>% 
    mt_add_pathways(in_col = "KEGG", out_col = "humancyc_db", pw_species = "hsapiens", pw_name = "humancyc", n_cpus = 5)
}


