#' Heatmap of Multiple Statistical Objects
#'
#' Generate a heatamp for a list of statistical results using \code{pheatmap::pheatmap}.
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_list List of stat names to plot; default NA (for all stat objects in D)
#' @param color_formula Expression used to color adjusted p-values; default: sign(fc) * -log10(p.adj)
#' @param met_anno rowData() column name to use as column annotation
#' @param signif_mark Expression used to indicate significant p-values
#' @param cluster_rows cluster by rows? default: F
#' @param cluster_cols cluster by cols? default: F
#' @param filter_signif show only metabolites significant in one or more results? default: T
#' @param show_colnames should column (metabolite) labels be included in the plot? default: F
#'
#' @return results$output: heatmap of statistical results
#'
#' @examples
#' \dontrun{D <- D %>% mt_plots_pheatmap_multstats(met_anno = "SUPER_PATHWAY", cluster_cols = T, cluster_rows = T)}
#'
#' @author KC
#'
#' @export

mt_plots_pheatmap_multstats <- function(D,
                                        stat_list=NA,
                                        color_formula=sign(statistic) * -log10(p.adj),
                                        met_anno,
                                        signif_mark,
                                        cluster_rows = F,
                                        cluster_cols = F,
                                        filter_signif = T,
                                        show_colnames = F){

  if (missing(signif_mark)) stop("signif_mark must be specified")
  color_formula <- dplyr::enquo(color_formula)
  signif_mark <- dplyr::enquo(signif_mark)

  pheat_arg <- list()
  pheat_arg$cluster_rows <- cluster_rows
  pheat_arg$cluster_cols <- cluster_cols
  pheat_arg$show_colnames <- show_colnames

  ## NTS: add checks for all parameters
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))


  # if NA, stat_list is all stat results
  if(is.na(stat_list)){
    stat_list <- MetaboTools:::mti_res_get_stats_entries(D)%>%
      purrr::map_chr(~.x$output$name)
  }

  if(length(stat_list) == 0){
    stop("There are no statistical results.")
  }

  # extract stat results
  stat_res <- stat_list %>% lapply(function(i){MetaboTools:::mti_get_stat_by_name(D, i) %>%
      dplyr::mutate(stat_name = i)
  })

  ## NTS: add a check to make sure the variables in the expression formulas are present

  # get stat names
  stat_names <- stat_res %>% sapply(function(x){
    stat_name = x$stat_name[1]
  }) %>% unname()

  # generate matrix to plot
  data_plot <- stat_res %>% lapply(function(x){
    stat_name = x$stat_name[1]
    x <- x %>%  dplyr::mutate(color = !!color_formula) %>%
      dplyr::mutate(stars = ifelse(!!signif_mark, "*", "")) %>%
      dplyr::select(., var, color, stars)
    colnames(x)[2] <- stat_name
    colnames(x)[3] <- paste0(stat_name, "_stars")
    x
  }) %>% purrr::reduce(dplyr::full_join, by="var")

  # filter significant metabolites?
  if(filter_signif == T){
    data_plot <- data_plot %>% dplyr::filter_at(dplyr::vars(dplyr::ends_with("_stars")), dplyr::any_vars(. == "*"))
  }

  # add annotations?
  if(!missing(met_anno)){
    rd <- rowData(D) %>% as.data.frame()
    tmp <- rd[,c("COMP_IDstr", met_anno)] %>% as.data.frame()
    colnames(tmp)[1] <- "var"
    tmp <- tmp %>% dplyr::semi_join(data_plot, by="var")
    rownames(tmp) <- tmp$var
    tmp <- tmp %>% dplyr::select(all_of(met_anno))

    pheat_arg$annotation_col <- tmp
  }

  # generate significance label table
  label_cols <- colnames(data_plot)[colnames(data_plot) %in% stat_names == F]
  data_label <- data_plot[,label_cols] %>% t()
  colnames(data_label) <- data_label[1,]
  data_label <- data_label[-1,] %>% as.data.frame()
  indx <- sapply(data_label, is.factor)
  data_label[indx] <- lapply(data_label[indx], function(x) as.character(x))
  rownames(data_label) <- gsub("_stars", "", rownames(data_label))
  pheat_arg$display_numbers <- data_label
  pheat_arg$fontsize_number <- 20

  # format data for plotting
  data_plot <- data_plot[,c("var", stat_names)] %>% t()
  colnames(data_plot) <- data_plot[1,]
  data_plot <- data_plot[-1,] %>% as.data.frame()
  indx <- sapply(data_plot, is.factor)
  data_plot[indx] <- lapply(data_plot[indx], function(x) as.numeric(as.character(x)))
  pheat_arg$mat <- data_plot

  # plot pheatmap
  re <- do.call(pheatmap::pheatmap, pheat_arg)

  # NTS: add option to convert to ggplot

  # add status information & plot
  funargs <- MetaboTools:::mti_funargs()
  metadata(D)$results %<>%
    MetaboTools:::mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Pheatmap of stat results"),
      output = list(re)
    )

  # return
  D

}
