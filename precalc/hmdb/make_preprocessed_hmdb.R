# MetaboTools external
#
# Parse HMDB xml file
#
# Parse the .xml file that can be downloaded from HMDB
#
# last update: 2019-01-14
# authors: Parviz Gomari
#

# dependencies
library(tidyverse)
library(XML)

# Download hmdb_metabolites.zip -------------------------------------------

hmdb_link <- "http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
zip_path <- str_split(hmdb_link, "\\/") %>% unlist() %>% tail(1)

download.file(hmdb_link, zip_path)
unzip(zip_path)
file.remove(zip_path)


# path to the saved hmdb xml file
xml_path <- zip_path %>% str_replace(".zip", ".xml")


# Preprocess parsed xml file and save as an .rds file ---------------------

# select the columns to retain for the preprocessing
retain_cols <- c(
  "accession", 
  "secondary_accessions", 
  "pathways"
  # "name", 
  # "iupac_name", 
  # "cas_registry_number",
  # "kegg_id"
  # "pubchem_compound_id", 
  # "chebi_id",
  # "chemspider_id", 
  # "bigg_id"
  # "biocyc_id",
  # "drugbank_id", 
  # "drugbank_metabolite_id", 
  # "foodb_id"
)

hmdb_raw <- xmlToDataFrame(xml_path)

hmdb_version <- hmdb_raw$version[[1]] %>% as.character()

hmdb_preprocessed <- 
  hmdb_raw %>% 
  dplyr::select(retain_cols) %>% 
  # convert all factors to characters
  mutate_if(is.factor, as.character) %>% 
  # remove empty pathways
  filter(pathways != "\n  ") %>% 
  rowwise() %>% 
  mutate(
    secondary_accessions = 
      str_replace_all(secondary_accessions, "HMDB", "|HMDB") %>% 
      str_remove("^\\|"),
    
    SMP = 
      str_match_all(pathways, "(SMP[0-9]{5})(map[0-9]{5})?") %>% 
      .[[1]] %>% 
      .[,1] %>% 
      list(),
    pathway_names = 
      str_split(pathways, "SMP[0-9]{5}|map[0-9]{5}") %>% 
      unlist() %>% 
      .[. != ""] %>% 
      list()
  ) %>% 
  unnest(SMP, pathway_names) %>% 
  # separate SMPDB pathway IDs from KEGG IDs
  mutate(
    KEGG = str_match(SMP, "map[0-9]{5}"),
    SMP = str_replace(SMP, "map[0-9]{5}", "")
  ) %>% 
  unnest(HMDB_id = strsplit(secondary_accessions, "\\|")) %>% 
  # remove entries with empty secondary accessions
  filter(HMDB_id != "\n  ") %>% 
  dplyr::select(-pathways) %>% 
  dplyr::select(
    accession, 
    HMDB_id, # Note: in the datasets, the secondary accessions are used
    SMP, 
    KEGG,
    pathway_name = pathway_names) %>% 
  distinct()


# Delete xml file and save preprocessed file
file.remove(xml_path)
saveRDS(hmdb_preprocessed, glue("hmdb/hmdb_preprocessed_{hmdb_version}.rds"))
