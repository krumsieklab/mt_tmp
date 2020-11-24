#' Wrapper function to run global statistics visualization
#'
#' Generating PCA, UMAP and heatmap plots of dataset
#'
#' @param D \code{SummarizedExperiment} input
#' @param annos_pca_umap List of names in colData column to be colored in PCA and UMAP, column names should be quoted by quote()
#' @param do_pca Logical. Whether to generate PCA plot
#' @param pca_options A list of parameters for \code{mt_plots_PCA}, expressions should be quoted by quote()
#' @param do_umap Logical. Whether to generate UMAP plot 
#' @param umap_options A list of parameters for \code{mt_plots_UMAP}, expressions should be quoted by quote()
#' @param do_heatmap Logical. Whether to generate heatmap
#' @param heatmap_options A list of parameters for \code{mt_plots_pheatmap}, expressions should be quoted by quote()
#'
#' @return D with PCA, UMAP and heatmap plots
#' @examples
#' \dontrun{... %>% mtw_global_stats(annos_pca_umap = list(quote(Cohort)),
#' pca_options = list(scaledata = T, size = 2, ggadd=quote(scale_size_identity())),
#' umap_options = list(scaledata = T, size = 2, ggadd=quote(scale_size_identity())),
#' heatmap_options = list(scaledata = T, fontsize = 5)) %>% ...}
#'
#' @author Zeyu Wang
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#' @export

mtw_global_stats <-
  function(D,
           annos_pca_umap=NULL,
           do_pca = T,
           pca_options = list(),
           do_umap = T,
           umap_options = list(),
           do_pheatmap = T,
           heatmap_options = list()) {
    stopifnot("SummarizedExperiment" %in% class(D))
    
    # function to merge default and user input
    map_lists <- function(def, user) {
      # check if there are any illegal entries
      bad_args <- setdiff(names(user), names(def))
      if (length(bad_args) > 0) {
        stop(sprintf(
          "Invalid argument(s): %s",
          paste0(bad_args, collapse = ", ")
        ))
      }
      # now fill up entries in def
      def[names(user)] <- user
      # return
      def
    }
    
    pca_options_def = list(
      title = "PCA",
      scaledata = F,
      PCa = 1,
      PCb = 2,
      show = "scores",
      labelby = "",
      textrepel = T,
      ellipse = NA,
      expvarplot = F,
      store.matrices = F,
      ggadd = NULL,
      size = 2
    )
    
    umap_options_def = list(
      title = "UMAP",
      scaledata = F,
      labelby = "",
      textrepel = T,
      store.matrices = F,
      ggadd = NULL,
      n_neighbors = 15,
      size = 2
    )
    
    heatmap_options_def = list(
      scaledata = T,
      sym0 = F,
      ggadd = NULL,
      return.gg = T,
      gg.scale = 1,
      fontsize = 5
    )
    
    
    # PCA
    if(do_pca == T){
      pca_options <- map_lists(pca_options_def, pca_options)
      if(!missing(annos_pca_umap)){
        for(s in annos_pca_umap){
          pca_options$D <- D
          pca_options$color <- s
          D <- do.call("mt_plots_PCA", pca_options)
        }
      }else{
        pca_options$D <- D
        D <- do.call("mt_plots_PCA", pca_options)
      }
    }
    
    # UMAP
    if(do_umap == T){
      umap_options <- map_lists(umap_options_def, umap_options)
      if(!missing(annos_pca_umap)){
        for(s in annos_pca_umap){
          umap_options$D <- D
          umap_options$color <- s
          D <- do.call("mt_plots_UMAP",umap_options)
        }
      }else{
        umap_options$D <- D
        D <- do.call("mt_plots_UMAP",umap_options)
      }
    }
    # Heatmap
    if(do_pheatmap == T){
      heatmap_options <- map_lists(heatmap_options_def,heatmap_options)
      heatmap_options$D <- D
      D <- do.call("mt_plots_pheatmap",heatmap_options)
    }
    D
  }
