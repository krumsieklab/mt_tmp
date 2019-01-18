# MetaboTools
#
# Add pathway information
#
# Adds custom pathways to the already existing SummarizedExperiment
# data structure using a flatfile.
#
# last update: 2019-01-17
# author: Parviz Gomari
#

# todo: document output

# dependencies
library(tidyverse)
library(data.table)
library(readxl)

# main
mt_anno_pathways_from_file <- function(
  D,                  # SummarizedExperiment input
  in_col,             # column to use for pathway fetching. The selected column must contain metabolite identifiers (e.g. HMBD, KEGG, ChEBI, etc)
  out_col,            # a new column name for D to output pathway information to
  flatfile,           # path where the pathway annotation flat file is stored
  sheet,              # sheet name or number or number to read in flat file
  met_ID_col,         # flatfile colname: this column should contain metabolite IDs
  pw_ID_col,          # flatfile colname: this column should contain pathway IDs
  pw_name_col,        # flatfile colname: this column should contain pathway names
  export_raw_db       # OPTIONAL. Export the pathway database to a directory. Must be a string containing the path name with a .xlsx extension.
) {
  
  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  if (missing(flatfile))
    stop("flatfile must be given to fetch pathway annotation file from")
  
  if (!file.exists(flatfile))
    stop(glue::glue("{flatfile} does not exist. input a valid flatfile path."))
  
  if(!in_col %in% names(rowData(D)))
    stop(glue::glue("in_col is invalid. please enter an existing column name."))
  
  
  mti_logstatus(glue::glue("reading annotation file: {basename(flatfile)}, sheet: {sheet}"))
  # check flatfile colnames
  pwdb <- read_excel(path=flatfile, sheet=sheet)
  
  flatfile_cols <- c(met_ID_col, pw_ID_col, pw_name_col)
  valid_names <- flatfile_cols %in% names(pwdb)
  if (!all(valid_names)) {
    invalid_names <- flatfile_cols[!valid_names]
    stop(sprintf("Non-existent flatfile column names. Please replace: %s \n with one of: %s", 
                 str_c(invalid_names, collapse = ", "),
                 str_c(names(pwdb), collapse = ", ")))
  }
  
  pwdb %<>% 
    dplyr::select(met_ID = !!met_ID_col, 
                  ID = !!pw_ID_col,
                  pathway_name = !!pw_name_col)
  
  # create a dataframe that enables the mapping between pathway
  # names and IDs. Included also are num_total, num_measured,
  # num_pw_total, and num_pw_measured (see below for further details)
  # in this dataframe map
  
  # num_total - overall number of metabolites in entire database (this will 
  # be a redundant, repeating number, identical in every row… but it’s the 
  # easiest way to store it right now)
  # pwdb$accession is used here, since there is a many-to-one mapping 
  # between HMDB_id (secondary accessions) and accession
  num_total <- 
    pwdb$met_ID %>% 
    unique() %>% 
    length()
  
  # num_measured - overall number of measured metabolites found in the DB
  num_measured = 
    pwdb$met_ID %>% 
    # for this intersect, similar to the one in pwdb_summary (below),
    # it is assumed that the IDs in the dataset D is non-redundant.
    # else this number may not be accurate.
    intersect(rowData(D)[[in_col]]) %>% 
    length()
  
  mti_logstatus(glue::glue("summarizing pathway information"))
  pwdb_summary <- pwdb
  # using methods from data.table reduces runtime by almost 10x as compared
  # to dplyr
  pwdb_summary <- 
    setDT(pwdb_summary)[, `:=`(
      
      num_total = num_total,
      num_measured = num_measured,
      
      # num_pw_total - the total number of metabolites in that pathway 
      # here, accession is used.
      # (overall DB background)
      num_pw_total = uniqueN(met_ID), 
      
      # num_pw_measured - the number of measured metabolites in that pathway 
      # (for the M type of analysis on the actual measured background).
      num_pw_measured = 
        sum(met_ID %in% rowData(D)[[in_col]], na.rm = TRUE)
    ),
    by = ID] %>% 
    # some IDs might have more than two names, however, these will be discarded
    # for now
    unique(by = c("ID")) %>% 
    subset(!is.na(ID), 
           select = -met_ID) 
  
  
  mti_logstatus(glue::glue("nesting pathway IDs"))
  # nest all the pathway IDs given our lieblings input IDs 
  pwdb_reduced <- 
    pwdb %>%
    group_by(met_ID) %>% 
    filter(!is.na(ID)) %>% 
    distinct(met_ID, ID) %>% 
    nest(ID, .key = ID) %>% 
    mutate(ID = 
             ID %>%
             unlist(recursive = FALSE) %>% 
             as.list())
  
  
  # match the nested pathways to our lieblings IDs
  match_idx <- 
    match(rowData(D)[[in_col]],
          pwdb_reduced$met_ID)
  
  pw_col <- pwdb_reduced$ID[match_idx]
  
  # add the pathway IDs into D
  rowData(D)[[out_col]] <- pw_col
  
  # add pathway map to the metadata of D
  metadata(D)$pathways[[out_col]] <- 
    pwdb_summary
  
  
  if(!missing(export_raw_db)) {
    openxlsx::write.xlsx(pwdb, export_raw_db)
  }
  
  
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('added pathway annotations using %s', basename(flatfile))
    )
  
  D
}



if (FALSE) {
  # Example -----------------------------------------------------------------
  mt_logging(console=T) 
  D_alone <- 
    mt_files_load_metabolon(codes.makepath("packages/metabotools/sampledata.xlsx"), "OrigScale") %>% 
    mt_anno_pathways_from_file(
      in_col = "HMDb_ID", out_col = "smp_db", 
      flatfile = codes.makepath("packages/metabotools_external/hmdb/hmdb_preprocessed_4.0.xlsx"),
      sheet = "hmdb",
      met_ID_col = "HMDB_id", pw_ID_col = "SMP", pw_name_col = "pathway_name"
    ) 
}

