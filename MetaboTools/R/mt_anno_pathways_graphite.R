#' Add pathway information.
#'
#' Adds custom pathways to the already existing SummarizedExperiment
#'  data structure using the graphite package.
#'
#' @param D \code{SummarizedExperiment} input
#' @param in_col rowData olumn to use for pathway fetching. The selected column must contain metabolite identifiers (e.g. HMBD, KEGG, ChEBI, etc)
#' @param out_col new column name for rowData to output pathway information to
#' @param pwdb_species name of the species the data was measured in. Use pathwayDatabases() after loading graphite to see a full list of databases and species
#' @param pwdb_name name of the pathway database to use. Use pathwayDatabases() after loading graphite to see a full list of databases and species
#' @param n_cores number of cores to use for parallelizaion (used in convertIdentifiers)
#' @param export_raw_db # OPTIONAL. Export the pathway database to a directory. Must be a string containing the path name with a .xlsx extension.
#'
#' @return rowData: new pathway annotation column added
#' @return $result$pathways: a dataframe of pathway information
#'
#' @examples
#' # annotate metabolites using the humancyc database from the graphite package
#' \dontrun{... %>%
#'     mt_add_pathways(in_col = "KEGG",
#'                     out_col = "humancyc_db",
#'                     pwdb_species = "hsapiens",
#'                     pwdb_name = "humancyc",
#'                     n_cores = 5) %>%
#' ...}
#'
#' @author Parviz Gomari
#'
#' @export

# main
mt_anno_pathways_graphite <- function(
  D,             # SummarizedExperiment input
  in_col,        # column to use for pathway fetching. The selected column must contain metabolite identifiers (e.g. KEGG, ChEBI, HMBD, etc)
  out_col,       # name of the column to output pathway information to in D
  pwdb_species,    # name of the species the data was measured in. Use pathwayDatabases() after loading graphite to see a full list of databases and species
  pwdb_name,       # the name of the pathway database to use. Use pathwayDatabases() after loading graphite to see a full list of databases and species
  n_cores = 2,    # number of cores to use for parallelizaion (used in convertIdentifiers)
  export_raw_db  # OPTIONAL. Export the pathway database to a directory. Must be a string containing the path name with a .xlsx extension.
) {

  # check arguments, SummarizedExperiment, and exactly one cutoff argument must be non-NA
  stopifnot("SummarizedExperiment" %in% class(D))


  ######################################################################################
  ## This will be relocated into another script where all ID combinations
  ## have been remapped to graphite databases.
  ######################################################################################
  # to see all available databasaes and species, use:
  # pathwayDatabases() # %>% filter(species == "hsapiens")
  pwdb <- graphite::pathways(species = pwdb_species, database = pwdb_name)

  # use given number of cores, this should be part of input
  options(Ncpus = n_cores)

  # convert from ChEBI IDs (given by graphite) to our labelings ID
  mti_logstatus(glue::glue("converting graphite ChEBI IDs to {in_col} IDs"))
  con_id <- dplyr::if_else(in_col == "KEGG", "KEGGCOMP", in_col)
  pwdb <- pwdb %>% graphite::convertIdentifiers(con_id)

  # graphite comes with as a list of pathways given a database
  # given this list, subselect the metabolite entries per pathway
  mti_logstatus(glue::glue("subselecting pathways that contain metabolites"))
  pwdb <-
    mcmapply(function(pwname) {

      pw <- pwdb[[pwname]]
      pw %>%
        graphite::edges(which = "metabolites") %>%
        dplyr::mutate(pathway_name = pwname,
               ID = graphite::pathwayId(pw))

    },
    names(pwdb),
    SIMPLIFY = F,
    mc.cleanup = T,
    mc.cores = n_cores)

  # convert the pathway list into a dataframe
  pwdb <-
    pwdb %>%
    dplyr::bind_rows() %>%
    dplyr::select(src, dest, pathway_name, ID) %>%
    tidyr::gather(key = tmp, value = mappingID, -c(pathway_name, ID)) %>%
    dplyr::select(mappingID, pathway_name, ID) %>%
    tidyr::drop_na() %>% # ensure that there are no NAs in any rows
    dplyr::distinct()

  # create a dataframe that enables the mapping between pathway
  # names and IDs. Included also are num_total, num_measured,
  # num_pw_total, and num_pw_measured (see below for further details)
  pwdb_summary <-
    pwdb %>%
    dplyr::mutate(
      # num_total - overall number of metabolites in entire database (this will
      # be a redundant, repeating number, identical in every row… but it’s the
      # easiest way to store it right now)
      num_total = dplyr::n_distinct(mappingID),
      # num_measured - overall number of measured metabolites found in the DB
      num_measured =
        mappingID %>%
        intersect(rowData(D)[[in_col]]) %>%
        length()
    ) %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(
      # num_pw_total - the total number of metabolites in that pathway
      # (overall DB background)
      num_pw_total = dplyr::n(),
      # num_pw_measured - the number of measured metabolites in that pathway
      # (for the M type of analysis on the actual measured background).
      num_pw_measured =
        mappingID %>%
        intersect(rowData(D)[[in_col]]) %>%
        length()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mappingID) %>%
    dplyr::distinct()

  # nest all the pathway IDs given our labelings input IDs
  pwdb_reduced <-
    pwdb %>%
    dplyr::group_by(mappingID) %>%
    dplyr::summarise(IDs = stringr::str_c(ID, collapse = ", ")) %>%
    dplyr::mutate(IDs = stringr::str_split(IDs, ", "))

  ######################################################################################
  ## End of move
  ######################################################################################

  # have to do this to take a variable input name for the columns
  # used for joining in the next step
  left_index <- dplyr::enquo(in_col)
  right_index <- "mappingID"
  by <-  purrr::set_names(dplyr::quo_name(right_index), dplyr::quo_name(left_index))

  # match the nested pathways to our labelings IDs
  pw_col <- D %>%
    rowData %>%
    .[in_col] %>%
    as.data.frame() %>%
    dplyr::left_join(pwdb_reduced, by = by) %>%
    .$IDs

  # add the pathway IDs into D
  rowData(D)[[out_col]] <- pw_col

  # add pathway map to the metadata of D
  metadata(D)$pathways[[out_col]] <- pwdb_summary


  if(!missing(export_raw_db)) {
    openxlsx::write.xlsx(pwdb, export_raw_db)
  }


  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('added pathway annotations using the %s database', pwdb_name)
    )

  D
}
