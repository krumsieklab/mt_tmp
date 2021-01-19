#' Creates KEGG pathways visualization
#'
#' ADD DESCRIPTION
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_name name of the statistics object to apply metab_filter to
#' @param gene_id string name of the rowData column containing the gene identifiers.
#' @param met_id String name of the rowData column containing the metabolite identifiers. If this parameter given, cpd.data cannot be used.
#' @param metab_filter if given, filter will be applied to data and only variables satisfying the condition will be included
#' @param color_scale if given, this will be used to map colors to a continuous scale
#' @param color_range numeric (positive), if given, indicates the color range (-color_range, +color_range). If missing, color_range will be determined internally.
#' @param show_only_filtered logical, if TRUE generate pathway list only based on filtered variables, otherwise pathways will be generated based on all variables.
#' @param n_pathways (optional) number of pathways to output. Most populated pathway will be plotted first.
#' @param path_database character, the directory path of KEGG pathway data file (.xml) and image file (.png). If the path does not exist, the function will create it. Default path_database = "./Pathview_database" (subfolder in the current working directory).
#' @param path_output character, the directory path of the function output files. If the path does not exist, the function will create it. Default path_output ="./Pathview_output" (subfolder in the current working directory).
#' @param add_pwname_suffix logical, if TRUE will add the pathway name to the output filename. If FALSE, will use what stored in out.suffix for all files. Default add.pwname.to.filename=FALSE.
#' @param db_path path to the pathway database to read annotations from.
#' @param gene.data either vector (single sample) or a matrix-like data (multiple sample). Vector should be numeric with gene IDs as names or it may also be character of gene IDs. Character vector is treated as discrete or count data. Matrix-like data structure has genes as rows and samples as columns. Row names should be gene IDs. Here gene ID is a generic concepts, including multiple types of gene, transcript and protein uniquely mappable to KEGG gene IDs. KEGG ortholog IDs are also treated as gene IDs as to handle metagenomic data. Check details for mappable ID types. Default gene.data=NULL.
#' @param cpd.data The same as gene.data, excpet named with IDs mappable to KEGG compound IDs. See \code{pathview::pathview} for details.
#' @param low,mid,high each is a list of two colors with "gene" and "cpd" as the names. This argument specifies the color spectra to code gene.data and cpd.data. Default spectra (low-mid-high) "green"-"gray"-"red" and "yellow"-"gray"-"blue" are used for gene.data and cpd.data respectively. The values for 'low, mid, high' can be given as color names ('red'), plot color index (2=red), and HTML-style RGB, ("\#FF0000"=red).
#' @param pathway.id character vector, the KEGG pathway ID(s), usually 5 digit, may also include the 3 letter KEGG species code. If missing, the function will find all KEGG pathway annotations for the given KEGG identifiers.
#' @param same.layer logical, controls if node colors are to be plotted in the same layer as the pathway graph. If FALSE, output generation will be faster, but output plots will be larger in size.
#' @param out.suffix character, the suffix to be added after the pathway name as part of the output graph file. Default out.suffix="pathview".
#' @param min.nnodes
#' @param expand.node
#' @param split.group
#' @param map.symbol
#' @param node.sum
#' @param \dots  see \code{pathview::pathview} for additional pathview arguments
#'
#' @return $result: pathview images
#'
#' @examples
#' \dontrun{# plot all pathways with at least one significant metabolite from the statistical comparison "comp" in them
#' mt_plots_pathview(D = D,
#'                   met_id="KEGG_mapped",
#'                   stat_name = "comp",
#'                   color_scale = -sign(fc)*log10(p.adj),
#'                   color_range = -log10(0.01),
#'                   metab_filter = p.adj < 0.05,
#'                   show_only_filtered = TRUE,
#'                   path_database = "./Pathview_database",
#'                   path_output = "./results/pathview",
#'                   same.layer = F,
#'                   add_pwname_suffix = T
#'                   ) %>%
#'                   ...}
#'
#' @author Elisa Benedetti
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import pathview
#' @import SummarizedExperiment
#'
#' @export
mt_plots_pathview <- function(D,
                              stat_name,
                             # only one between gene_id and gene.data can be given
                             gene_id = NULL,
                             gene.data = NULL,
                             # only one between met_id and cpd.data can be given
                             met_id = NULL,
                             cpd.data = NULL,
                             # if gene_id or met_id is given, variables can be selected from the results of a statistical analysis
                             metab_filter,
                             color_scale,
                             color_range,
                             # only show pathways that include filtered results?
                             show_only_filtered = FALSE,
                             # colors for genes and metabolite data
                             low = list(gene = "green", cpd = "yellow"),
                             mid =list(gene = "gray", cpd = "gray"),
                             high = list(gene = "red", cpd ="blue"),
                             # pathways to plot
                             pathway.id,
                             n_pathways,
                             # result directories
                             path_database = "./Pathview_database",
                             path_output = "./Pathview_output",
                             # set to false for speed-up (output files will be bigger in size though)
                             same.layer = TRUE,
                             out.suffix = "pathview",
                             add_pwname_suffix = FALSE,
                             db_path = system.file("extdata", "precalc/pathview/KeggPathways.Rds", package = "MetaboTools"),
                             # other pathview::pathview arguments
                             min.nnodes = 1,
                             expand.node = TRUE,
                             split.group = TRUE,
                             map.symbol = FALSE,
                             node.sum = "median"
                             ) {

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  ## if both gene_id and gene.data are given, throw an error
  if(!is.null(gene_id) & !is.null(gene.data))
    stop("The function can only use either gene.data or gene_id. Please provide only one of the two.")
  ## if both met_id and cpd.data are given, throw an error
  if(!is.null(met_id) & !is.null(cpd.data))
    stop("The function can only use either cpd.data or met_id. Please provide only one of the two.")
  ## if gene_id is provided, check that is is a valid column name of the rowData
  if(!is.null(gene_id)) {
    if (!(gene_id %in% colnames(rowData(D))))
      stop(sprintf("%s is not contained in the rowData of the Summarized Experiment", gene_id))}
  ## if met_id is provided, check that is is a valid column name of the rowData
  if(!is.null(met_id) ) {
     if(!(met_id %in% colnames(rowData(D))))
      stop(sprintf("%s is not contained in the rowData of the Summarized Experiment", met_id))}
  ## if gene_id is provided, check that it is length 1
  if(!is.null(gene_id)) {
    if(length(gene_id)!=1)
      stop(sprintf("%s can only be a single column", gene_id))}
  ## if met_id is provided, check that it is length 1
  if(!is.null(met_id)) {
    if(length(met_id)!=1)
      stop(sprintf("%s can only be a single column", met_id))}
  ## if metab_filter is given, stat_name must also be given and either met_id or gene_id must be given as well
  if(!missing(metab_filter)) {
    if(missing(stat_name))
      stop("In order to use metab_filter, stat_name must be given")
    if(is.null(gene_id) & is.null(met_id))
      stop("In order to use metab_filter, one betweeen gene_id and met_id must be given")}
  ## if n.pathway is given, it must be numeric
  if(!missing(n_pathways)) {
    if(class(n_pathways)!="numeric")
      stop("n_pathways must be numeric")}
  ## if show_only_filtered is TRUE, metab_filter must be given
  if(show_only_filtered) {
    if(missing(metab_filter))
      stop("show_only_filtered can be TRUE only if metab_filter is given")}
  ## in order for add_pwname_suffix to work when TRUE, pathway.id must be missing
  if(add_pwname_suffix){
    if(!(missing(pathway.id)))
      stop("add_pwname_suffix can only be TRUE if pathway.id is missing")
  }

  ## rowData
  rd <- rowData(D) %>%
    as.data.frame() %>%
    dplyr::mutate(var = rownames(D))

  ## stat
  if(!missing(stat_name)){
    stat <- mti_get_stat_by_name(D, stat_name) %>%
      dplyr::inner_join(rd, by = "var")
  }else{
    stat <- rd
  }

  ## Add color variable
  if(!missing(color_scale)){
    color_scale_q <- dplyr::enquo(color_scale)
    # add color variable according to input
    stat <- stat %>%
      dplyr::mutate(color=!!color_scale_q)
  } else {
    # if not given, set color to 1
    stat <- stat %>%
      dplyr::mutate(color=1)
  }

  ## FILTER METABOLITES
  if(!missing(metab_filter)){
    metab_filter_q <- dplyr::enquo(metab_filter)
    # filter results
    var <- stat %>%
      dplyr::filter(!!metab_filter_q)
    # collect variable names of filtered results
    var <- var$var
    # set color variable of filtered out variables to 0
    stat$color[!(stat$var %in% var)] <- 0
  }

  # if gene_id is provided, extract identifiers from the rowData
  if(!is.null(gene_id)) {
    # remove duplicated identifiers if they occur
    stat <- stat[!duplicated(stat[[gene_id]]),]
    # create pathview variable
    gene.data <- data.frame(color=stat$color[!is.na(stat[[gene_id]])])
    # remove rows with NAs in the identifiers
    rownames(gene.data) <- stat[[gene_id]][!is.na(stat[[gene_id]])]
  }
  # if met_id is provided, extract identifiers from the rowData
  if(!is.null(met_id)) {
    # remove duplicated identifiers if they occur
    stat <- stat[!duplicated(stat[[met_id]]),]
    # create pathview variable
    cpd.data <- data.frame(color=stat$color[!is.na(stat[[met_id]])])
    # remove rows with NAs in the identifiers
    rownames(cpd.data) <- stat[[met_id]][!is.na(stat[[met_id]])]
  }

  # if path_output is provided, check if the folder exists, otherwise create it
  if (!file.exists(path_output)){
    dir.create(path_output)
  }
  # if path_database is provided, check if the folder exists, otherwise create it
  if (!file.exists(path_database)){
    dir.create(path_database)
  }

  # if no pathway.id list is provided, find annotations for kegg identifiers
  if(missing(pathway.id)) {

    # load KEGG pathway database
    load(db_path)

    # build one big dataframe with all pathway informations
    pwdf <- do.call(rbind, pwdb)

    if (!is.null(gene.data)) {
      if(!is.null(rownames(gene.data))) {
        ids <- rownames(gene.data)
        if(show_only_filtered) {
          ids <- ids[ids %in% rownames(gene.data)[gene.data$color != 0]]
        }
      } else {
        ids <- gene.data
      }

      if(length(ids) != 0) {
        # find gene pathway annotations
        g_anno <- lapply(ids, function(x) {
          pwdf$ID[pwdf$src==x] %>% unique()
        })
        names(g_anno) <- ids
        # build one long list
        g_anno_list <- do.call(c, g_anno)
        # find most common pathway for genes
        pw_gene <- g_anno_list %>% table() %>% as.data.frame()
        colnames(pw_gene) <- c("pathway","Freq")
        # pathway list ordered according to the number of genes with that annotation
        pw_gene <- pw_gene[order(pw_gene$Freq,decreasing = TRUE),]
        # find names of these pathways
        g_pw_names <- lapply(pw_gene$pathway, function(x) {
          pwdf$name[pwdf$ID==x] %>% unique()
        })
        names(g_pw_names) <- pw_gene$pathway
        # build one long list
        pw_names <- do.call(c, g_pw_names)
        # remove ":" from pathway ids for pathview
        pw_gene$pathway <- gsub(":", "", pw_gene$pathway)
        names(pw_names) <- gsub(":", "", names(pw_names))
        pw <- pw_gene
      } else {
        warning("Filtering returned an empty matrix")
        pw <- list()
        pw$pathway <- NULL
        gene.data <- NULL
      }
      if(!(is.null(pw$pathway))) {
        # save pathway list only if list of variables to output is not empty
        pathway.id <- pw$pathway
        }
    }

    if (!is.null(cpd.data)) {
      if(!is.null(rownames(cpd.data))) {
        ids <- rownames(cpd.data)
        if(show_only_filtered) {
          ids <- ids[ids %in% rownames(cpd.data)[cpd.data$color != 0]]
        }
      } else {
        ids <- cpd.data
      }

      if(length(ids)!=0) {
        # find metabolite pathway annotations
        m_anno <- lapply(ids, function(x) {
          pwdf$ID[pwdf$dest==x] %>% unique()
          # cbind(pwdf$ID[pwdf$dest==x] %>% unique(),pwdf$name[pwdf$dest==x] %>% unique())
        })
        names(m_anno) <- ids
        # build one long list
        m_anno_list <- do.call(c, m_anno)
        # find most common pathway for metabolites
        pw_met <- m_anno_list %>% table() %>% as.data.frame()
        colnames(pw_met) <- c("pathway","Freq")
        # pathway list ordered according to the number of metabolites with that annotation
        pw_met <- pw_met[order(pw_met$Freq,decreasing = TRUE),]
        # find names of these pathways
        m_pw_names <- lapply(pw_met$pathway, function(x) {
          pwdf$name[pwdf$ID==x] %>% unique()
        })
        names(m_pw_names) <- pw_met$pathway
        # build one long list
        pw_names <- do.call(c, m_pw_names)
        # remove ":" from pathway ids for pathview
        pw_met$pathway <- gsub(":", "", pw_met$pathway)
        names(pw_names) <- gsub(":", "", names(pw_names))
        pw <- pw_met
      } else {
        warning("Filtering returned an empty matrix")
        pw <- list()
        pw$pathway <- NULL
        cpd.data <- NULL
      }
      if(!(is.null(pw$pathway))) {
        # save pathway list only if list of variables to output is not empty
        pathway.id <- pw$pathway
        }
    }

    if (!is.null(gene.data) & !is.null(cpd.data)) {
      # find most common pathway for both genes and metabolites
      pw_list <- c(m_anno_list,g_anno_list)
      pw <- pw_list %>% table() %>% as.data.frame()
      colnames(pw) <- c("pathway","Freq")
      # pathway list ordered according to the number of metabolites/genes with that annotation
      pw <- pw[order(pw$Freq,decreasing = TRUE),]
      # find names of these pathways
      pw_names <- lapply(pw$pathway, function(x) {
        pwdf$name[pwdf$ID==x] %>% unique()
      })
      names(pw_names) <- pw$pathway
      # build one long list
      pw_names <- do.call(c, pw_names)
      # remove ":" from pathway ids for pathview
      pw$pathway <- gsub(":", "", pw$pathway)
      names(pw_names) <- gsub(":", "", names(pw_names))

      if(!(is.null(pw$pathway))) {
        # save pathway list only if list of variables to output is not empty
        pathway.id <- pw$pathway
        }
      }
  }

  ## Set color scale limits
  if(!missing(color_range)) {
    limit = list(gene=color_range, cpd=color_range)
  } else {
    # if(!is.null(met_id)) {
    #   limit$cpd = max(ceiling(abs(stat$color)), na.rm = TRUE)
    # }
    # if(!is.null(gene_id)) {
    #   limit$gene = max(ceiling(abs(stat$color)), na.rm = TRUE)
    # }
    # limit=list(gene=max(ceiling(abs(stat$color))), cpd=max(ceiling(abs(stat$color))))
    if(!is.null(cpd.data)) {
      if(!is.null(rownames(cpd.data))) {
        limit$cpd = max(ceiling(abs(cpd.data)), na.rm = TRUE)
      }
    }
    if(!is.null(gene.data)) {
      if(!is.null(rownames(gene.data))) {
        limit$gene = max(ceiling(abs(gene.data)), na.rm = TRUE)
      }
    }
  }
  if(limit$gene == 0) {limit$gene = 1}
  if(limit$cpd == 0) {limit$cpd = 1}

  # move working directory to kegg.dir (otherwise some files will be saved in the working directory even if another directory is provided)
  wd <- getwd()
  setwd(path_database)
  save.path <- getwd()
  setwd(wd)
  setwd(path_output)

  if(missing(pathway.id)) {
    file.create(paste0(getwd(),"/NO_RESULTS_AFTER_FILTERING.txt",sep=""))
    pv.out <- NULL
  } else {

    # removing problematic pathways from list if present -> they throw a weird error
    pathway.id <- pathway.id[!(pathway.id %in% c("mmu05206","mmu04666","mmu04723"))]

    if(!missing(n_pathways)) {
      if(n_pathways>length(pathway.id))
        warning(sprintf("n.pathway is %i, but there are only %i pathways, so %i pathways will be used", n_pathways, length(pathway.id), length(pathway.id)))
      pathway.id <- pathway.id[1:min(n_pathways,length(pathway.id))]
      pw_names <- pw_names[1:min(n_pathways,length(pathway.id))]
    }

    # print(sprintf("%d pathways detected", length(pathway.id)))

    lapply(1:length(pathway.id), function(x) {
      # print(sprintf("pathway %d, %s", x, pathway.id[x]))
      suppressMessages(
        pv.out <- pathview::pathview(gene.data = gene.data, cpd.data = cpd.data, pathway.id = pathway.id[x], kegg.dir = save.path,
                                     min.nnodes = min.nnodes, expand.node = expand.node, split.group = split.group, map.symbol = map.symbol,
                                     map.cpdname = map.cpdname, node.sum = node.sum,
                                     limit = limit, bins = bins,
                                     both.dirs = both.dirs, trans.fun = trans.fun, low = low, mid = mid, high = high, na.col = na.col,
                                     same.layer = same.layer, out.suffix = out.suffix)
      )
      # add pathway rank and name to filename
      if(add_pwname_suffix) {
        fname <- list.files(".",pattern=pathway.id[x])
        # isolate pathway name from filename
        m <- pw_names[names(pw_names)==substr(fname[length(fname)], 1, 8)]
        m <- gsub('[[:punct:]]+','',m)
        m <- stringr::str_replace_all(m," ","_")
        file.rename(from=fname,to=sub(pattern=sprintf("%s.%s",pathway.id[x], out.suffix),replacement=sprintf("%d_%s.%s",x,m,pathway.id[x]),fname))
      }
    }) %>% invisible()
  }

  setwd(wd)

  n_pw <- ifelse((!missing(pathway.id)),length(pathway.id),0)
  nn <- ifelse((!missing(pathway.id)),pw_names,NA)

  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("%d Pathway images saved in %s", n_pw, path_output),
      output = NULL,
      output2 = nn
    )

  # return
  D
}
