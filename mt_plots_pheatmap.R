#' 
#' Heatmap plot by pheatmap::pheatmap 
#' 
#' 
#' @param D summarized experiment object
#' @param scaledata scaling the data, TRUE by default
#' @param sym0 make color scale symmetric around 0? (should only be used for scaled data), default: F 
#' @param fD function to transform/scale \code{t(assay(D))}, ie \code{mat = fD(t(assay(D)))} will be plotted

#' @param return.gg should pheatmap object be converted to gg object, TRUE for default.
#' @param gg.scale scaling of plot to be converted to gg object
#' @param ggadd  further elements/functions to add (+) to the ggplot object

#' @param \dots  see \code{pheatmap::pheatmap} for pheatmap arguments 
#' @return object \code{SummarizedExperiment}, see \code{metabotools} conventions for the details
#' @note all \code{pheatmap::pheatmap} arguments can be passed 
#' @author mubu
#' @references \code{\link{https://github.com/raivokolde/pheatmap}}
#' @keywords ~heatmap ~pheatmap
#' @examples 
#' 
#' D %>%
#' mt_plots_pheatmap(annotation_row = c("SUPER_PATHWAY", "PLATFORM", "RI"), 
#'                   annotation_col = c("GROUP_DESC","BATCH_MOCK","gender"), 
#'                   fD = function(x) scale(exp(scale(x))),
#'                   clustering_distance_cols =  "correlation",
#'                   clustering_distance_rows = "minkowski")
#' 
#' @export


mt_plots_pheatmap <- function(D, scaledata=F, sym0=F, fD = function(x){ if(scaledata) return(scale(x)); x}, # metabotools arguments
                              
                              # pheatmap::pheatmap arguments
                              color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, 
                              border_color = "grey60", cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, 
                              clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", 
                              clustering_callback = pheatmap:::identity2, cutree_rows = NA, cutree_cols = NA,  
                              treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows, 50, 0), 
                              treeheight_col = ifelse((class(cluster_cols) == "hclust") || cluster_cols, 50, 0), 
                              legend = TRUE, legend_breaks = NA, legend_labels = NA, 
                              annotation_row = NA, annotation_col = NA, annotation_colors = NA, annotation_legend = TRUE, 
                              annotation_names_row = TRUE, annotation_names_col = TRUE, drop_levels = TRUE, show_rownames = T, show_colnames = T, 
                              main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = "%.2f", 
                              number_color = "grey30", fontsize_number = 0.8 * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL, 
                              labels_col = NULL, filename = NA, width = NA, height = NA, silent = TRUE, na_col = "#DDDDDD",
                              ggadd=NULL,  
                              # returned plot type, and specs
                              return.gg = T, gg.scale = 1, gg.ymin = 1 - gg.scale, gg.xmin = 1 - gg.scale, 
                              gg.xmax = gg.scale, gg.ymax = gg.scale,  
                              ...){
  
  # upon Jan's comment annotation_col and annotation_row are swapped for compatibility with SummarizedExperiment
  
  # get all inputs 
  aa = c(as.list(environment()), list(...))
  
  # fD(t(assay(D))) will be heatmapped
  x = t(assay(D))
  # if x doesn't have rownames, use numbers 1:nrow
  if (is.null(rownames(x))) rownames(x) <- 1:nrow(x)
  # store
  x.colnames = colnames(x)
  x.rownames = rownames(x)
  aa$mat = fD(x)
  
  # keep only pheatmap::pheatmap parameters
  aa = aa[!(names(aa) %in% c("D","scaledata","fD","return.gg", "gg.scale", "gg.ymin", "gg.xmin", "gg.xmax", "gg.ymax"))]
  # annotations will be added later
  aa = aa[!(names(aa) %in% c("annotation_col", "annotation_row"))]
  
  # deprecated pheatmap parameter 'annotation', see pheatmap::pheatmap
  aa$annotation = NA
  
  # annotate rows with given variables in attr(D, "colData")
  if(!is.na(annotation_col[1])){
    annotation_col = D %>% colData %>% as.data.frame %>% `[`(,annotation_col,drop = F)
    rownames(annotation_col) = x.rownames
    aa$annotation_row = annotation_col
  }
  
  # annotate columns with given variables in attr(D, "rowData")
  if(!is.na(annotation_row[1])){
    annotation_row = D %>% rowData %>% as.data.frame %>% `[`(,annotation_row,drop = F)
    rownames(annotation_row) = x.colnames
    aa$annotation_col = annotation_row
  }
  
  # symmetric around zero?
  if (sym0) {
    cap <- max(abs(aa$mat)) # from -cap to +cap
    cs = length(color) # number of color steps
    aa$breaks = seq(-cap, cap, 2*cap/cs) # technical stuff
  }
  
  # plot pheatmap 
  re <- do.call(pheatmap::pheatmap, aa)
  
  # cast pheatmap object to gg object 
  if(return.gg){
    re <- ggplot(data.frame(x = 0:1, y = 0:1), aes_(x = ~x, y = ~y)) +
      geom_blank() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
      annotation_custom(re$gtable, xmin = gg.xmin, xmax = gg.xmax, ymin = gg.ymin, ymax = gg.ymax) +
      theme_void()
    
    # add custom elements?
    if (!is.null(ggadd)) re <- re+ggadd
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


