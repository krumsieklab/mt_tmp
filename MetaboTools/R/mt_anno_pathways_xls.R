#' Add pathway information
#'
#' Adds custom pathways to the already existing SummarizedExperiment data structure using a flat file.
#' NOTE: The flat file should be an Excel file, containing 3 columns:
#' \itemize{
#'   \item{met_id} contains metabolite IDs
#'   \item{pw_id} contains pathway IDs
#'   \item{pw_name} contains pathway names}
#' The pathway columns pw_id and pw_name should be in long format, meaning if a metabolite is assigned to multiple pathways,
#' each metabolite pathway pair should appear in a separate row.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param in_col rowData column to use for pathway fetching. The selected column must contain metabolite
#'    identifiers (e.g. HMBD, KEGG, ChEBI, etc).
#' @param out_col New column name for rowData to output pathway information to.
#' @param file Path where the pathway annotation flat file is stored.
#' @param sheet Sheet name or number to read in flat file.
#' @param met_id Name of flat file column containing metabolite IDs.
#' @param pw_id Name of flat file colname containing pathway IDs.
#' @param pw_name Name of flat file colname containing pathway names.
#' @param raw_db_outfile OPTIONAL. Name of file to export the pathway database to. Must be a string containing the path name with a
#'    .xlsx extension.
#'
#' @return rowData: New pathway annotation column added.
#' @return $results$pathways: A dataframe of pathway information.
#'
#' @examples
#' \dontrun{# annotate metabolites using the SMP column of the pathway database flat file
#' ... %>%
#'       mt_anno_pathways_from_file(in_col = "HMDb_ID",
#'                                  out_col = "janpw",
#'                                  file = codes.makepath("packages/metabotools_external/hmdb/hmdb_preprocessed_4.0.csv"),
#'                                  met_id = "HMDB_id",
#'                                  pw_id = "SMP",
#'                                  pw_name = "pathway_name") %>%
#' ...}
#'
#' @author PG
#'
#' @importFrom data.table :=
#'
#' @export
mt_anno_pathways_xls <- function(D,
                                 in_col,
                                 out_col,
                                 file,
                                 sheet,
                                 met_id,
                                 pw_id,
                                 pw_name,
                                 raw_db_outfile) {

  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if (missing(file))
    stop("file must be given to fetch pathway annotation file from")

  if (!file.exists(file))
    stop(glue::glue("{file} does not exist. input a valid flat file path."))

  if(!in_col %in% names(rowData(D)))
    stop(glue::glue("in_col is invalid. please enter an existing column name."))


  mti_logstatus(glue::glue("reading annotation file: {basename(file)}, sheet: {sheet}"))
  # check file colnames
  pwdb <- readxl::read_excel(path=file, sheet=sheet)

  flatfile_cols <- c(met_id, pw_id, pw_name)
  valid_names <- flatfile_cols %in% names(pwdb)
  if (!all(valid_names)) {
    invalid_names <- flatfile_cols[!valid_names]
    stop(sprintf("Non-existent flat file column names. Please replace: %s \n with one of: %s",
                 stringr::str_c(invalid_names, collapse = ", "),
                 stringr::str_c(names(pwdb), collapse = ", ")))
  }

  pwdb %<>%
    dplyr::select(met_id = !!met_id,
                  ID = !!pw_id,
                  pathway_name = !!pw_name)

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
    pwdb$met_id %>%
    unique() %>%
    length()

  # num_measured - overall number of measured metabolites found in the DB
  num_measured =
    pwdb$met_id %>%
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
    data.table::setDT(pwdb_summary)[, `:=`(

      num_total = num_total,
      num_measured = num_measured,

      # num_pw_total - the total number of metabolites in that pathway
      # here, accession is used.
      # (overall DB background)
      num_pw_total = data.table::uniqueN(met_id),

      # num_pw_measured - the number of measured metabolites in that pathway
      # (for the M type of analysis on the actual measured background).
      num_pw_measured =
        sum(met_id %in% rowData(D)[[in_col]], na.rm = TRUE)
    ),
    by = ID] %>%
    # some IDs might have more than two names, however, these will be discarded
    # for now
    unique(by = c("ID")) %>%
    subset(!is.na(ID),
           select = -met_id)


  mti_logstatus(glue::glue("nesting pathway IDs"))
  # nest all the pathway IDs given our lieblings input IDs
  pwdb_reduced <-
    pwdb %>%
    dplyr::group_by(met_id) %>%
    dplyr::filter(!is.na(ID)) %>%
    dplyr::distinct(met_id, ID) %>%
    tidyr::nest(ID, .key = ID) %>%
    dplyr::mutate(ID =
             ID %>%
             unlist(recursive = FALSE) %>%
             as.list())


  # match the nested pathways to our lieblings IDs
  match_idx <-
    match(rowData(D)[[in_col]],
          pwdb_reduced$met_id)

  pw_col <- pwdb_reduced$ID[match_idx]

  # add the pathway IDs into D
  rowData(D)[[out_col]] <- pw_col

  # add pathway map to the metadata of D
  metadata(D)$pathways[[out_col]] <-
    pwdb_summary


  if(!missing(raw_db_outfile)) {
    openxlsx::write.xlsx(pwdb, raw_db_outfile)
  }


  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('added pathway annotations using %s', basename(file))
    )

  D
}
