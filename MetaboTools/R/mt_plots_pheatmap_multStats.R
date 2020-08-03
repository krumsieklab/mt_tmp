#' Heatmap of Multiple Statistical Objects
#'
#' Generate a heatamp for a list of statistical results using \code{pheatmap::pheatmap}.
#'
#' @param D \code{SummarizedExperiment} input
#' @param stat_list List of stat names to plot; default NA (for all stat objects in D)
#' @param color_formula Expression used to color adjusted p-values; default: sign(fc) * -log10(p.adj)
#' @param color_signif Boolean. Color only significant values? default: F
#' @param color_sigval a vector of length 2 containing positive and negative thresholds of significance
#' @param met_anno rowData() column name to use as column annotation
#' @param signif_expr Expression used to indicate significant p-values
#' @param signif_mark Expression used to indicate significant p-values
#' @param cluster_rows cluster by rows? default: F
#' @param cluster_cols cluster by cols? default: F
#' @param filter_signif show only metabolites significant in one or more results? default: T
#' @param show_colnames should column (metabolite) labels be included in the plot? default: F
#' @param color_signif only color significant values? default: F
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
                                        signif_expr,
                                        signif_mark = "•",
                                        mark_size = 20,
                                        show_mark = T,
                                        color_signif = F,
                                        color_sigval,
                                        filter_signif = F,
                                        cluster_rows = F,
                                        cluster_cols = F,
                                        show_colnames = F,
                                        main = "Stat Results Overview"){

  color_formula <- dplyr::enquo(color_formula)

  pheat_arg <- list()
  pheat_arg$cluster_rows <- cluster_rows
  pheat_arg$cluster_cols <- cluster_cols
  pheat_arg$show_colnames <- show_colnames
  pheat_arg$main <- main

  ## NTS: add checks for all parameters
  # validate and extract arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  if(color_signif==T){
    if(missing(color_sigval)){
      stop("color_signif is set to True, but no values provided for color_sigval")
    }else if(length(color_sigval) != 2){
      stop("color_sigval must consist of two and only two elements")
    }else if(!all(is.numeric(color_sigval))){
      stop("color_sigval must contain only numeric values")
    }
  }

  if(missing(signif_expr)){
    if(filter_signif == T){
      stop("filter_signif is TRUE, but no value for signif_expr was provided.")
    }else if(color_signif == T){
      stop("color_signif is TRUE, but no value for signif_expr was provided.")
    }
  }else{
    signif_expr <- dplyr::enquo(signif_expr)
  }


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

  if(!missing(signif_expr)){
    # determine which metabolites are significant if signif_expr provided

    # generate matrix to plot
    data_plot <- stat_res %>% lapply(function(x){
      stat_name = x$stat_name[1]
      x <- x %>%  dplyr::mutate(color = !!color_formula) %>%
        dplyr::mutate(signif_col = ifelse(!!signif_expr, signif_mark, "")) %>%
        dplyr::select(., var, color, signif_col)
      colnames(x)[2] <- stat_name
      colnames(x)[3] <- paste0(stat_name, "_signif_col")
      x
    }) %>% purrr::reduce(dplyr::full_join, by="var")

    all_plot <- data_plot

    # filter significant metabolites
    if(filter_signif == T){
      data_plot <- data_plot %>% dplyr::filter_at(dplyr::vars(dplyr::ends_with("_signif_col")), dplyr::any_vars(. == signif_mark))
    }

    # generate significance label table
    if(show_mark==T){
      label_cols <- colnames(data_plot)[colnames(data_plot) %in% stat_names == F]
      data_label <- data_plot[,label_cols] %>% t()
      colnames(data_label) <- data_label[1,]
      data_label <- data_label[-1,] %>% as.data.frame()
      indx <- sapply(data_label, is.factor)
      data_label[indx] <- lapply(data_label[indx], function(x) as.character(x))
      rownames(data_label) <- gsub("_signif_col", "", rownames(data_label))
      pheat_arg$display_numbers <- data_label
      pheat_arg$fontsize_number <- mark_size
    }

  }else{
    data_plot <- stat_res %>% lapply(function(x){
      stat_name = x$stat_name[1]
      x <- x %>%  dplyr::mutate(color = !!color_formula) %>%
        dplyr::select(., var, color)
      colnames(x)[2] <- stat_name
      x
    }) %>% purrr::reduce(dplyr::full_join, by="var")
  }

  # add annotations?
  if(!missing(met_anno)){
    rd <- rowData(D) %>% as.data.frame()
    tmp <- rd[,c("name", met_anno)] %>% as.data.frame()
    colnames(tmp)[1] <- "var"
    tmp$var <- make.names(tmp$var)
    tmp <- tmp %>% dplyr::semi_join(data_plot, by="var")
    rownames(tmp) <- tmp$var
    tmp <- tmp %>% dplyr::select(all_of(met_anno))

    pheat_arg$annotation_col <- tmp
  }

  # format data for plotting
  data_plot <- data_plot[,c("var", stat_names)] %>% t()
  colnames(data_plot) <- data_plot[1,]
  data_plot <- data_plot[-1,] %>% as.data.frame()
  indx <- sapply(data_plot, is.factor)
  data_plot[indx] <- lapply(data_plot[indx], function(x) as.numeric(as.character(x)))
  pheat_arg$mat <- data_plot

  # if coloring by significance, need to calculate breaks
  if(color_signif == T){

    # generate colors and color breaks
    pheat_arg$color <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                       "RdYlBu")))(100)
    n = length(pheat_arg$color)
    color_breaks <- seq(min(pheat_arg$mat, na.rm=T), max(pheat_arg$mat, na.rm=T), length.out = n + 1)

    # find breaks that represent significant values
    left_thresh <- min(color_sigval)
    right_thresh <- max(color_sigval)

    flat_mat <- sort(unlist(pheat_arg$mat))
    left_break <- max(flat_mat[flat_mat < left_thresh])
    right_break <- min(flat_mat[flat_mat > right_thresh])

    sig_left_idx <- color_breaks < left_thresh
    sig_right_idx <- color_breaks > right_thresh
    sig_left_idx[max(which(sig_left_idx))+1] <- TRUE

    sig_idx <- Reduce("|", list(sig_left_idx, sig_right_idx))

    # change threshold breaks so no significant values are caught in between
    color_breaks[max(which(sig_left_idx))+1] <- left_break
    color_breaks[min(which(sig_right_idx))] <- right_break

    pheat_arg$color[sig_idx[-c(1)]==F] <- "#D3D3D3"

    pheat_arg$breaks <- color_breaks

  }


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
