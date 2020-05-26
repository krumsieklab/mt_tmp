#' Create KEGG identifiers from HMDB identifiers.
#' 
#' Creates annotations and merges them into current SummarizedExperiment. 
#' Performs "left-joins", i.e. leaves the original SE unchanged and just adds information where it can be mapped.
#'
#' @param D \code{SummarizedExperiment} input
#' @param col_input string name of the rowData column containing the HMDB identifiers
#' @param col_output string name of the new rowData column to be created with the KEGG identifiers
#'
#' @return rowData: new annotations added
#'
#' @examples
#' # Load data, sheet with sample annotations, sheet with metabolite annotations and add KEGG identifiers
#' D <- 
#'   # load raw data
#'   mt_files_data_xls(file=file, sheet="data", samplesInRows=T, ID="SAMPLE_NAME") %>% 
#'   # sample annotations from metabolomics run
#'   mt_files_anno_xls(file=file, sheet="sampleinfo", annosfor="samples", IDanno = "SAMPLE_NAME") %>% 
#'   # metabolite annotations
#'   mt_files_anno_xls(file=file, sheet="metinfo", annosfor="metabolites", IDanno="BIOCHEMICAL", IDdata = "name") %>% 
#'   # add KEGG identifiers
#'   mt_anno_metabolites_HMDB2KEGG(col_input="HMDb_id", col_output="KEGG") %>%
#'   ...
#' 
#' @author EB
#' 
#' @export

mt_anno_metabolites_HMDB2KEGG <- function(D,
                                          col_input,
                                          col_output = "KEGG"
                                          ) {
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(col_input))
    stop("col_input must be given")
  # check that col_input is a string
  if(!(class(col_input)=="character"))
    stop("col_input must be a string")
  # check that col_input is a valid column names of the rowData
  if(!(col_input %in% colnames(rowData(D))))
    stop(sprintf("%s is not contained in the rowData of the Summarized Experiment", col_input))
  if(col_output %in% colnames(rowData(D)))
    stop(sprintf("%s already exists, please choose another name", col_output))

  # load look-up table of HMDB to KEGG identifiers
  load(data.makepath("MT_precalc/pathview/MetaboliteMapping.Rds"))
  # get rowData
  df <- as.data.frame(rowData(D))
  # rename col_input for left_join use
  names(df)[which(names(df)==col_input)] <- "MappingIDs"
  
  # join with automatic metabolite mapping
  newdf <- df %>% 
    dplyr::left_join(MetaboliteMapping[,c("secondary_accessions","kegg_id")], by = c("MappingIDs" = "secondary_accessions")) %>%
    dplyr::left_join(MetaboliteMapping[,c("accession","kegg_id")], by=c("MappingIDs" = "accession")) %>% 
    distinct()
  # merge KEGG identifiers
  newdf[[col_output]] <- newdf$kegg_id.x
  newdf[[col_output]] [is.na(newdf[[col_output]])] <- newdf$kegg_id.y[is.na(newdf[[col_output]])]
  # drop redundant columns
  newdf <- subset(newdf, select=-c(kegg_id.x,kegg_id.y))
  
  # go back to original name
  # rename col_input for left_join use
  names(newdf)[which(names(newdf)=="MappingIDs")] <- col_input
  
  # overwrite entries stored in the manually curated file
  dd <- read.csv2(file=data.makepath("MT_precalc/pathview/MetaboliteMapping_manual.csv"), sep=",")
  dd$Name <- as.character(dd$Name)
  dd$HMDB <- as.character(dd$HMDB)
  dd$KEGG <- as.character(dd$KEGG)
  
  # find identifiers to overwrite
  ow <- intersect(newdf[[col_input]],dd$HMDB)
  # subselect
  dd <- dd[dd$HMDB %in% ow,]
  # get order of identifiers in rowData
  dd <- dd[match(newdf[[col_input]][which(newdf[[col_input]] %in% ow)], dd$HMDB),]
  # overwrite identifiers
  newdf[[col_output]][which(newdf[[col_input]] %in% ow)] <- dd$KEGG
  # update rowData
  rowData(D) <- DataFrame(newdf)
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded KEGG annotations for %i out of %i metabolites", length(newdf[[col_output]][!is.na(newdf[[col_output]])]), length(newdf[[col_output]]))
    )
  
  # return
  D
  
} 