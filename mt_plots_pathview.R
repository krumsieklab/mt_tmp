require(pathview)

#' mt_plots_pathview
#'
#' Creates KEGG pathways visualization
#'
#' @param D \code{SummarizedExperiment} input
#' @param gene.id string name of the rowData column containing the gene identifiers.
#' @param met.id string name of the rowData column containing the metabolite identifiers
#' @param statname name of the statistics object to apply metab_filter to
#' @param metab_filter if given, filter will be applied to data and only variables satisfying the condition will be included
#' @param gene.data either vector (single sample) or a matrix-like data (multiple sample). Vector should be numeric with gene IDs as names or it may also be character of gene IDs. Character vector is treated as discrete or count data. Matrix-like data structure has genes as rows and samples as columns. Row names should be gene IDs. Here gene ID is a generic concepts, including multiple types of gene, transcript and protein uniquely mappable to KEGG gene IDs. KEGG ortholog IDs are also treated as gene IDs as to handle metagenomic data. Check details for mappable ID types. Default gene.data=NULL.
#' @param cpd.data the same as gene.data, excpet named with IDs mappable to KEGG compound IDs. Over 20 types of IDs included in CHEMBL database can be used here. Check details for mappable ID types. Default cpd.data=NULL. Note that gene.data and cpd.data can't be NULL simultaneously.
#' @param pathway.id character vector, the KEGG pathway ID(s), usually 5 digit, may also include the 3 letter KEGG species code.
#' @param kegg.dir character, the directory path of KEGG pathway data file (.xml) and image file (.png). If the path does not exist, the function will create it. Default kegg.dir="." (current working directory).
#' @param \dots  see \code{pathview::pathview} for pathview arguments
#' @return $result: pathview images
#' 
#' @examples
#' # map given metabolite KEGG identifiers and given genes Entrez identifiers to the pathway "hsa00010" and save images in the in the current directory
#' ... %>%
#' mt_plots_pathview(cpd.data=c("C02787","C08521","C01043","C11496","C07111"), 
#'                   gene.data=c("100008586" "100008587" "100008588" "100008589" "100009613"),
#'                   pathway.id="hsa00010",
#'                   kegg.dir = ".") %>%
#' ...
#' 
#' # map metabolite KEGG identifiers contained in the rowData column called "KEGG.ID" to the pathways "hsa00010","hsa00020" and save images in a folder called "Pathview" in the current directory
#' ... %>%
#' mt_plots_pathview(met.id="KEGG.ID",
#'                   pathway.id=c("hsa00010","hsa00020"),
#'                   kegg.dir = "./Pathview") %>%
#' ...
#'
#' @author Elisa Benedetti
#' 
mt_plots_pathview <- function(D,
                             # only one between gene.id and gene.data can be given
                             gene.id = NULL,
                             gene.data = NULL,
                             # only one between met.id and cpd.data can be given
                             met.id = NULL,
                             cpd.data = NULL,
                             # if gene.id or met.id is given, variables can be selected from the results of a statistical analysis
                             statname,
                             metab_filter,
                             pathway.id,
                             kegg.dir = ".",

                             # other pathview::pathview arguments
                             species = "hsa", cpd.idtype = "kegg", gene.idtype = "entrez", gene.annotpkg = NULL, min.nnodes = 3, kegg.native = TRUE,
                             map.null = TRUE, expand.node = FALSE, split.group = FALSE, map.symbol = TRUE, map.cpdname = TRUE, node.sum = "sum", 
                             discrete=list(gene=FALSE,cpd=FALSE), limit = list(gene = 1, cpd = 1), bins = list(gene = 10, cpd = 10), 
                             both.dirs = list(gene = T, cpd = T), trans.fun = list(gene = NULL, cpd = NULL), low = list(gene = "green", cpd = "blue"), 
                             mid =list(gene = "gray", cpd = "gray"), high = list(gene = "red", cpd ="yellow"), na.col = "transparent"
                             ) {
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  ## if both gene.id and gene.data are given, throw an error
  if(!is.null(gene.id) & !is.null(gene.data))
    stop("The function can only use either gene.data or gene.id. Please provide only one of the two.")
  ## if both met.id and cpd.data are given, throw an error
  if(!is.null(met.id) & !is.null(cpd.data))
    stop("The function can only use either cpd.data or met.id. Please provide only one of the two.")
  ## if gene.id is provided, check that is is a valid column name of the rowData
  if(!is.null(gene.id)) { 
    if (!(gene.id %in% colnames(rowData(D))))
      stop(sprintf("%s is not contained in the rowData of the Summarized Experiment", gene.id))}
  ## if met.id is provided, check that is is a valid column name of the rowData
  if(!is.null(met.id) ) {
     if(!(met.id %in% colnames(rowData(D))))
      stop(sprintf("%s is not contained in the rowData of the Summarized Experiment", met.id))}
  ## if gene.id is provided, check that it is length 1
  if(!is.null(gene.id)) { 
    if(length(gene.id)!=1)
      stop(sprintf("%s can only be a single column", gene.id))}
  ## if met.id is provided, check that it is length 1
  if(!is.null(met.id)) {
    if(length(met.id)!=1)
      stop(sprintf("%s can only be a single column", met.id))}
  ## if metab_filter is given, statname must also be given and either met.id or gene.id must be given as well
  if(!missing(metab_filter)) {
    if(missing(statname))
      stop("In order to use metab_filter, statname must be given")
    if(is.null(gene.id) & is.null(met.id))
      stop("In order to use metab_filter, one betweeen gene.id and met.id must be given")}

  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    mutate(var = rownames(D))
  
  ## stat
  if(!missing(statname)){
    stat <- mti_get_stat_by_name(D, statname) %>%
      inner_join(rd, by = "var")
  }else{
    stat <- rd
  }
  
  ## FILTER METABOLITES
  if(!missing(metab_filter)){
    metab_filter_q <- enquo(metab_filter)
    stat <- stat %>%
      filter(!!metab_filter_q)
    mti_logstatus(glue::glue("filter metabolites: {metab_filter_q} [{nrow(stat)} remaining]"))
  }
  
  # if gene.id is provided, extract identifiers from the rowData
  if(!is.null(gene.id)) {
    gene.data <- stat[[gene.id]]
  }
  # if met.id is provided, extract identifiers from the rowData
  if(!is.null(met.id)) {
    cpd.data <- stat[[met.id]]
  }

  # if kegg.dir is provided, check if the folder exists, otherwise create it
  if (!file.exists(kegg.dir)){
    dir.create(kegg.dir)
  } 
  
  # move working directory to kegg.dir (otherwise some files will be saved in the working directory even if another directory is provided)
  wd <- getwd()
  setwd(kegg.dir)
  save.path <- getwd()
  
  suppressMessages(
  pv.out <- pathview::pathview(gene.data = gene.data, cpd.data = cpd.data, pathway.id = pathway.id, kegg.dir = ".",
                               species = species, cpd.idtype = cpd.idtype, gene.idtype = gene.idtype, gene.annotpkg = gene.annotpkg, min.nnodes = min.nnodes, kegg.native = kegg.native,
                               map.null = map.null, expand.node = expand.node, split.group = split.group, map.symbol = map.symbol, map.cpdname = map.cpdname, node.sum = node.sum, 
                               discrete = discrete, limit = limit, bins = bins, 
                               both.dirs = both.dirs, trans.fun = trans.fun, high = high, na.col = na.col)
  )
  setwd(wd)
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Pathway images saved in %s", save.path),
      output = NULL,
      output2 = pv.out
    )
  
  # return
  D
}