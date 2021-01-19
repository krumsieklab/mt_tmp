#' Heatmap plot by pheatmap::pheatmap
#'
#' all \code{pheatmap::pheatmap} arguments can be passed
#' \code{\href{https://github.com/raivokolde/pheatmap}{https://github.com/raivokolde/pheatmap}}
#'
#' @param D summarized experiment object
#' @param scale_data scaling the data, TRUE by default
#' @param sym_zero make color scale symmetric around 0? (should only be used for scaled data), default: F
#' @param fun function to transform/scale \code{t(assay(D))}, ie \code{mat = fun(t(assay(D)))} will be plotted
#' @param silent Don't draw the table? A pheatmap argument, MetaboTools uses a different default. Default: T.
#' @param return_gg should pheatmap object be converted to gg object, TRUE for default.
#' @param gg_scale scaling of plot to be converted to gg object
#' @param gg_ymin
#' @param gg_xmin
#' @param gg_xmax
#' @param gg_ymax
#' @param ggadd  further elements/functions to add (+) to the ggplot object
#' @param \dots  see \code{pheatmap::pheatmap} for pheatmap arguments
#'
#'  @return object \code{SummarizedExperiment}, see \code{metabotools} conventions for the details
#'
#' @author mubu
#'
#' @examples
#' \dontrun{D %>%
#' mt_plots_heatmap(annotation_row = c("SUPER_PATHWAY", "PLATFORM", "RI"),
#'                   annotation_col = c("GROUP_DESC","BATCH_MOCK","gender"),
#'                   fun = function(x) scale(exp(scale(x))),
#'                   clustering_distance_cols =  "correlation",
#'                   clustering_distance_rows = "minkowski")
#'  }
#'
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @import pheatmap
#' @import SummarizedExperiment
#'
#' @export

mt_plots_heatmap <- function(D,
                              scale_data=F,
                              sym_zero=F,
                              fun = function(x){ if(scale_data) return(scale(x)); x},
                              silent = TRUE,
                              ggadd=NULL,
                              return_gg = T,
                              gg_scale = 1,
                              gg_ymin = 1 - gg_scale,
                              gg_xmin = 1 - gg_scale,
                              gg_xmax = gg_scale,
                              gg_ymax = gg_scale,
                              ...){

  # upon Jan's comment annotation_col and annotation_row are swapped for compatibility with SummarizedExperiment

  # get all inputs
  aa = c(as.list(environment()), list(...))

  # fun(t(assay(D))) will be heatmapped
  x = t(assay(D))
  if (any(is.na(x))) stop("Data matrix for heatmap cannot contain NAs")

  # if x doesn't have rownames, use numbers 1:nrow
  if (is.null(rownames(x))) rownames(x) <- 1:nrow(x)
  # store
  x.colnames = colnames(x)
  x.rownames = rownames(x)
  aa$mat = fun(x)

  # keep only pheatmap::pheatmap parameters
  aa = aa[!(names(aa) %in% c("D","scale_data","fun","return_gg", "gg_scale", "gg_ymin", "gg_xmin", "gg_xmax", "gg_ymax"))]
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
  if (sym_zero) {
    cap <- max(abs(aa$mat)) # from -cap to +cap
    cs = length(color) # number of color steps
    aa$breaks = seq(-cap, cap, 2*cap/cs) # technical stuff
  }

  # plot pheatmap
  re <- do.call(pheatmap::pheatmap, aa)

  # cast pheatmap object to gg object
  if(return_gg){
    re <- ggplot(data.frame(x = 0:1, y = 0:1), aes_(x = ~x, y = ~y)) +
      geom_blank() + scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
      annotation_custom(re$gtable, xmin = gg_xmin, xmax = gg_xmax, ymin = gg_ymin, ymax = gg_ymax) +
      theme_void()

    # add custom elements?
    if (!is.null(ggadd)) re <- re+ggadd
  }

  # fix ggplot environment
  if (D %>% mti_get_setting("ggplot_fix")) re <- mti_fix_ggplot_env(re)


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


