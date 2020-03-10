library(tidyverse)
library(parallel)
library(xml2)

# data source: http://www.hmdb.ca/system/downloads/current/saliva_metabolites.zip
met_xml <- xml2::read_xml(codes.makepath("snippets/packages/metabotools_external/pathview/hmdb_metabolites.xml"))

# check namespace used by parsed document
# the namespace is required for all xpaths, or else they won't work
doc_ns <- xml_ns(met_xml) %>% names()
search_tag <- "metabolite"

parsed_xml <- met_xml %>% xml_find_all(glue::glue(".//{doc_ns}:{search_tag}"))

df_list <- 
  mclapply(1:length(parsed_xml),
           function(idx) {
             print(idx)
             tibble(
               accession = parsed_xml[idx] %>% xml_find_first(".//d1:accession") %>% xml_text(),
               secondary_accessions = parsed_xml[idx] %>% xml_find_all(".//d1:secondary_accessions//d1:accession") %>% xml_text(),
               name = parsed_xml[idx] %>% xml_find_first(".//d1:name") %>% xml_text(),
               
               # cas_registry_number = parsed_xml[idx] %>% xml_find_all(".//d1:cas_registry_number") %>% xml_text() %>% list(),
               # pubchem_compound_id = parsed_xml[idx] %>% xml_find_all(".//d1:pubchem_compound_id") %>% xml_text() %>% list(),
               # chemspider_id = parsed_xml[idx] %>% xml_find_all(".//d1:chemspider_id") %>% xml_text() %>% list(),
               kegg_id = parsed_xml[idx] %>% xml_find_all(".//d1:kegg_id") %>% xml_text() %>% list() #,
               # meta_cyc_id = parsed_xml[idx] %>% xml_find_all(".//d1:meta_cyc_id") %>% xml_text() %>% list(),
               # pdb_id = parsed_xml[idx] %>% xml_find_all(".//d1:pdb_id") %>% xml_text() %>% list(),
               # chebi_id = parsed_xml[idx] %>% xml_find_all(".//d1:chebi_id") %>% xml_text() %>% list(),
               # 
               # disease_name = parsed_xml[idx] %>% xml_find_all(".//d1:disease//d1:name") %>% xml_text(),
               # disease_reference = parsed_xml[idx] %>% xml_find_all(".//d1:disease//d1:references") %>% xml_text()
             )
           },
           mc.cores = 4,
           mc.cleanup = TRUE)


MetaboliteMapping <- 
  bind_rows(df_list) %>% 
  # unnest(cas_registry_number, .drop = FALSE) %>% 
  # unnest(pubchem_compound_id, .drop = FALSE) %>% 
  # unnest(chemspider_id, .drop = FALSE) %>%
  unnest(kegg_id, .drop = FALSE) # %>%
  # unnest(meta_cyc_id, .drop = FALSE) %>%
  # unnest(pdb_id, .drop = FALSE) %>%
  # unnest(chebi_id, .drop = FALSE)

MetaboliteMapping %>% tail(30)

save(KeggHmdb_MetMapping, codes.makepath("snippets/packages/mtabotools_external/pathview/MetaboliteMapping.rds"))

