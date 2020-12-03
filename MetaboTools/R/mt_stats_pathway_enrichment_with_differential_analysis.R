# WHY IS stat_name MISSING?
#' Pathway enrichment using different methods
#'
#' This pathway enrichment analysis script calculats differential metabolites either
#' by a classical t-test or the GAGE package and then performs pathway enrichment
#' by using Fisher's exact test or GAGE.
#'
#' Implemented approaches:
#' 1. GAGE: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-161
#' 2. Fisher's exact test. Will NOT scale() data before. Data matrix can have NAs.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param pw_col Column containing pathways IDs.
#' @param grp_col Column name of column containing grouping variables.
#' @param ctrl_grp Name(s) of control group(s) found in grp_col column.
#' @param case_grp Name(s) of case/phenotype group(s) found in grp_col column.
#' @param method One of: "gage", "fishers".
#'
#' @return $pathways$enrichment_results: a dataframe containing the pathway enrichment results
#'
#' @examples
#' \dontrun{%>% mt_stats_pathway_enrichment("kegg_db",
#'                                 grp_col = "Group",
#'                                 ctrl_grp = "Vehicle",
#'                                 case_grp = c("treatment2", "treatment1") %>%
#' }
#'
#'
#' @author PG
#'
#' @export
mt_stats_pathway_enrichment_with_differential_analysis <- function(D,
                                                                   pw_col,
                                                                   grp_col,
                                                                   ctrl_grp,
                                                                   case_grp,
                                                                   method = "gage") {

  stopifnot("SummarizedExperiment" %in% class(D))
  if (!method %in% c("fishers","gage")) stop("'method' must be either 'gage' or 'fishers'")

  meta_D <- metadata(D)

  if(!"pathways" %in% names(meta_D)) stop("'pathways' does not exist in current SummarizedExperiment input")

  # Check if given pathway column actually exists
  if (!pw_col %in% names(meta_D$pathways)) stop(sprintf("'%s' not found in metabolite annotations.", pw))

  pw_id_map <-
    meta_D$pathways[[pw_col]] %>%
    dplyr::distinct(ID, pathway_name)

  row.data <- rowData(D) %>% as.data.frame()
  D.df <- assay(D) %>% as.data.frame()
  if(!("COMP_IDstr" %in% names(row.data))) {
    # If there isnt compound id (eg. if data from WCM core), make our own
    row.data$COMP_IDstr <- paste("cid", seq(nrow(row.data)), sep="")

    D.df$name <- rownames(D.df)
    nr <- nrow(D.df)


    row.data$new.name <- rownames(row.data)
    D.df <- dplyr::inner_join(D.df, row.data %>% dplyr::select(new.name, COMP_IDstr),
                       by=c("name"="new.name"))
    row.data <- row.data %>% dplyr::select(-new.name)
    if(nr != nrow(D.df)) {
      stop("Number of rows changed after merge, something is wrong")
    }
    rownames(D.df) <- D.df$COMP_IDstr
    D.df <- D.df %>% dplyr::select(-c(COMP_IDstr, name))
  }


  geneset <-
    row.data %>%
    dplyr::select(COMP_IDstr, !!rlang::sym(pw_col)) %>%
    dplyr::filter(!!rlang::sym(pw_col) != "NULL") %>%
    tidyr::unnest(!!rlang::sym(pw_col)) %>%
    dplyr::filter(!!rlang::sym(pw_col) != "NA") %>%
    dplyr::distinct()

  pw_names <-
    geneset[[pw_col]] %>%
    unique()

  geneset_list <-
    sapply(
      pw_names,
      simplify = FALSE,
      USE.NAMES = TRUE,
      function(pw_name) {
        geneset %>%
          dplyr::filter(!!rlang::sym(pw_col) == pw_name) %>%
          .$COMP_IDstr
      }
    )

  # create indices for control and case groups
  # NOTE: it is assumed here that the order of sample groups in colData
  # follows the order of samples in the assay dataframe
  sample_indices <-
    colData(D) %>%
    as.data.frame() %>%
    dplyr::transmute(!!rlang::sym(grp_col),
                     sample_idx = dplyr::row_number())

  ctrl_samples <-
    sample_indices %>%
    dplyr::filter(!!rlang::sym(grp_col) %in% ctrl_grp) %>%
    .$sample_idx

  case_samples <-
    sample_indices %>%
    dplyr::filter(!!rlang::sym(grp_col) %in% case_grp) %>%
    .$sample_idx


  if (method == "gage") {

    # perform pathway enrichment using the GAGE package
    enrichment_results <-
      D.df %>%
      gage::gage(gsets = geneset_list,
           ref = ctrl_samples,
           samp = case_samples,
           same.dir = FALSE, # suggsted for KEGG pathways
           compare = "unpaired", # unpaired samples
           set.size = c(5, 500)) %>%
      .$greater %>%
      as.data.frame() %>%
      tibble::rownames_to_column("pathway_ID") %>%
      dplyr::select(pathway_ID,
                    p_value = p.val,
                    p_value_adjusted = q.val,
                    mean_foldchange = stat.mean) %>%
      dplyr::filter(!is.na(p_value)) %>%

      dplyr::left_join(pw_id_map, by = c("pathway_ID" = "ID")) %>%
      dplyr::select(pathway_name, dplyr::everything())


  } else if (method == "fishers") {

    # Algorithm summary:
    # - calculate metabolite p-values per group
    # - adjust p-values using FDR
    # - assign significance to adjusted p-value at 0.05 level
    # - perform Fisher's exact test
    # - adjust Fisher's exact test p-value
    # - calculate mean fold change based on mean log value of cases and ctrls

    met_stats <-
      D.df %>%
      t() %>%
      as.data.frame() %>%

      # assign sample index
      # NOTE: the assumption here is that each row is ordered
      # as the rows in colData
      dplyr::mutate(sample_idx = dplyr::row_number()) %>%

      # select samples
      dplyr::filter(sample_idx %in% c(ctrl_samples, case_samples)) %>%
      dplyr::mutate(grp = dplyr::if_else(sample_idx %in% ctrl_samples, "ctrl", "case")) %>%
      dplyr::select(-sample_idx) %>%

      # Calculate p-values and the group means
      dplyr::summarise_if(is.numeric, list(p_value = ~ stats::t.test(. ~ .data$grp)$p.value,
                                    mean_ctrl = ~ stats::t.test(. ~ .data$grp)$estimate %>% as.numeric() %>% .[1],
                                    mean_case = ~ stats::t.test(. ~ .data$grp)$estimate %>% as.numeric() %>% .[2])) %>%
      tidyr::gather(met_id, vals) %>%
      tidyr::separate(met_id, c("met_id", "val_type"), "_", extra = "merge") %>%
      tidyr::spread(val_type, vals) %>%
      dplyr::mutate(p_value_adjusted = stats::p.adjust(p_value, method = "fdr"))

    enrichment_results <-
      met_stats %>%

      # assign significance
      dplyr::mutate(significant = dplyr::if_else(p_value_adjusted < 0.05, TRUE, FALSE),
                    n_total = dplyr::n(),
                    n_total_sig = sum(significant)) %>%
      dplyr::inner_join(geneset, by = c("met_id" = "COMP_IDstr")) %>%
      dplyr::group_by(!!rlang::sym(pw_col)) %>%

      # calculate summary numbers for Fisher's test
      dplyr::summarise(n_total = unique(n_total),
                       n_total_sig = unique(n_total_sig),
                       n_pw = dplyr::n(),
                       n_pw_sig = sum(significant),
                       mean_case = mean(mean_case),
                       mean_ctrl = mean(mean_ctrl)) %>%
      dplyr::filter(n_pw >= 5) %>%

      # calculate contingency table entries
      dplyr::mutate(s_p = n_pw_sig,
                    ns_p = n_pw - n_pw_sig,
                    s_np = n_total_sig - n_pw_sig,
                    ns_np = n_total - (s_p + ns_p + s_np)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(p_value =
                      stats::fisher.test(matrix(c(s_p, s_np, ns_p, ns_np), nrow = 2)) %>%
                      .$p.value) %>%
      dplyr::ungroup() %>%
      dplyr::rename(ID = !!rlang::sym(pw_col)) %>%
      dplyr::left_join(pw_id_map, by = "ID") %>%
      dplyr::transmute(pathway_name,
                       pathway_ID = ID,
                       p_value,
                       p_value_adjusted = stats::p.adjust(p_value, method = "fdr"),
                       mean_foldchange = mean_case - mean_ctrl) %>%
      dplyr::arrange(p_value)

  } else {

    stop("Bug.")

  }


  metadata(D)$pathways$enrichment_results <-
    tibble::as_tibble(enrichment_results)


  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('performed pathway enrichment on %s pathways using "%s"',
                       nrow(enrichment_results), method)
    )

  D
}
