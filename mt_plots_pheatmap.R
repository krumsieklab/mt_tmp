library(pheatmap)

#' all pheatmap inputs identical to phetamap::pheatmap function 
#' except for 
#' D : summarized experiment object
#' mat : if mat is not given, mat = assay(D)
#' annotation_by_SummarizedExperiment: if True, 
#'         columns are annotated with elementMetaData in D with given column names annotation_col
#'         rows are annotated with colData in D with given column names annotation_row
#' return.gg: if TRUE, pheatmep object will be cast to gg object and returned, otherwise, pheatmsp object returns 
#' gg.scale: scaling of the plot while casting to gg object 
#' gg.(ymin,ymax,xmin,xmax): min,max of gg coordinates if gg.return = T
#' for all other input argumets see pheatmap::pheatmap
#' 
mt_plots_pheatmap <- function(D, scaledata=F,  mat = NULL, color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, 
                             border_color = "grey60", cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, 
                             clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", 
                             clustering_callback = pheatmap:::identity2, cutree_rows = NA, cutree_cols = NA,  
                             treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows, 50, 0), 
                             treeheight_col = ifelse((class(cluster_cols) == "hclust") || cluster_cols, 50, 0), 
                             legend = TRUE, legend_breaks = NA, legend_labels = NA, 
                             annotation_row = NA, annotation_col = NA, annotation = NA, annotation_colors = NA, annotation_legend = TRUE, 
                             annotation_names_row = TRUE, annotation_names_col = TRUE, drop_levels = TRUE, show_rownames = T, show_colnames = T, 
                             main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = "%.2f", 
                             number_color = "grey30", fontsize_number = 0.8 * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL, 
                             labels_col = NULL, filename = NA, width = NA, height = NA, silent = TRUE, na_col = "#DDDDDD",
                             return.gg = T, gg.scale = 1, gg.ymin = 1 - gg.scale, gg.xmin = 1 - gg.scale, gg.xmax = gg.scale, gg.ymax = gg.scale, 
                             annotation_by_SummarizedExperiment= T, ...){
  
  
  # get all inputs 
  aa = c(as.list(environment()), list(...))
  
  # if no external datagiven to be plotted, assay(D) will be heatmapped
  if(is.null(mat)){
    aa$mat = t(assay(D))
    aa = aa[-1]  
  }
  
  # scale?
  if (scaledata) aa$mat <- scale(aa$mat)
  
  # annotate columns with given variables in attr(D, "elementMetadata")
  if(annotation_by_SummarizedExperiment & !is.na(annotation_col[1])){
    annotation_col = attr(D,"elementMetadata") %>% as.data.frame %>% `[`(,annotation_col,drop = F)
    rownames(annotation_col) = colnames(aa[[1]])
    aa$annotation_col = annotation_col
  }
  
  # deprecated pheatmap parameter 'annotation', see pheatmap::pheatmap  
  # annotate columns with given variables in attr(D, "elementMetadata")
  if(annotation_by_SummarizedExperiment & !is.na(annotation[1])){
    annotation = attr(D,"elementMetadata") %>% as.data.frame %>% `[`(,annotation,drop = F)
    rownames(annotation) = colnames(aa[[1]])
    aa$annotation = annotation
  }
  
  # annotate rowss with given variables in attr(D, "colData")
  if(annotation_by_SummarizedExperiment & !is.na(annotation_row[1])){
    annotation_row = attr(D,"colData") %>% as.data.frame %>% `[`(,annotation_row,drop = F)
    rownames(annotation_row) = rownames(aa[[1]])
    aa$annotation_row = annotation_row
  }
  
  # plot pheatmap 
  re <- do.call(pheatmap, aa)
  
  # cast pheatmap object to gg object 
  if(return.gg){
    re <- ggplot(data.frame(x = 0:1, y = 0:1), aes_(x = ~x, y = ~y)) +
      geom_blank() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
      annotation_custom(re$gtable, xmin = gg.xmin, xmax = gg.xmax, ymin = gg.ymin, ymax = gg.ymax) +
      theme_void()
  }
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Pheatmap..."),
      output = list(re)
    )
  
  # return
  D
}



if(FALSE){
  mt_logging(console=T) 
  D <- 
    mt_files_load_metabolon(codes.makepath("packages/metabotools/sampledata.xlsx"), "OrigScale") %>%
    mt_plots_PCA_mult(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>%
    mt_plots_sampleboxplot() %>%
    mt_plots_qc_missingness() %>%
    mt_pre_filtermiss(metMax=0.2) %>%
    mt_pre_filtermiss(sampleMax=0.1) %>%
    mt_pre_batch_median(batches = "BATCH_MOCK") %>%
    mt_pre_norm_quot() %>%
    mt_pre_trans_log() %>%
    mt_plots_qc_dilutionplot(comp="num1") %>%
    mt_plots_qc_dilutionplot(comp="Group") %>%
    mt_pre_impute_knn() %>%
    mt_plots_sampleboxplot(color=Group) %>%
    mt_plots_PCA(color=Group, shape=BATCH_MOCK, size=NUM_MOCK) %>% 
    mt_plots_pheatmap(annotation_col = c("SUPER_PATHWAY", "PLATFORM", "RI"), 
                      annotation_row = c("GROUP_DESC","BATCH_MOCK","gender"), return.gg = T)
  
  metadata(D)$results[[15]]$output
  
}
