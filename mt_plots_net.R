require(GGally)
require(sna)
require(visNetwork)

#' Network plotting
#'
#' @param D \code{SummarizedExperiment} input
#' @param statname name of the test to take correlations from
#' @param corr_filter filter for correlation values to plot
#' @param node_coloring name of the test to use for node coloring
#' 
#' @return assay: not altered
#' @return $result: interactive network plot
#'
#' @examples
#' #' # in the context of a SE pipeline
#' ... %>% mt_plots_net(statsname = "xxx") %>% ...    # standard call
#' ... %>% mt_plots_net(statsname = "xxx", corr_filter = p.adj < 0.5, node_coloring="Li's") %>% ...    # filters only significant correlations and colors the nodes according to the results in the indicated test
#'
#' @author EB
#' @export
#'

mt_plots_net <- function(
  D,                               # SummarizedExperiment input
  statname,                        # name of the correlation matrix to plot
  corr_filter = p.value < 0.05,    # filter
  node_coloring                    # name of the statistical test to use for node coloring
){
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(statname))
    stop("statname must be given to plot the network")
    
  ## rowData
  rd1 <- subset(rowData(D), select=which(names(rowData(D))=="name")) %>%
        as.data.frame() %>%
        mutate(var1 = rownames(D))
  colnames(rd1)[colnames(rd1)=="name"] <- "name1"
  rd2 <- subset(rowData(D), select=which(names(rowData(D))=="name")) %>%
    as.data.frame() %>%
    mutate(var2 = rownames(D))
  colnames(rd2)[colnames(rd2)=="name"] <- "name2"
  
  ## stat
  data_plot <- mti_get_stat_by_name(D, statname) %>%
    inner_join(rd1, by = "var1") %>%
    inner_join(rd2, by = "var2")
  
  ## define node attributes
  ids <- unique(data_plot[,which(colnames(data_plot)=="name1")])
  labels <- ids
  
  if(!(missing(node_coloring))) {
    test <- mti_get_stat_by_name(D, node_coloring);
    test <- test[match(ids, test$var),];
    map <- sign(test$statistic)*log10(test$p.value);
    node_color <- colorRampPalette(c('blue', 'white', 'red'))(length(map))[rank(map)]
  } else node_color <- "darkblue"

  nodes <- data.frame(id=ids, label=labels, color=node_color)
  
  ## apply filter on correlations
  if(!missing(corr_filter)){
    mti_logstatus("filter correlations")
    corr_filter_q <- enquo(corr_filter)
    data_plot <- data_plot %>%
      filter(!!corr_filter_q)
  }
  
  ## define edge attributes
  # rescale correlation to [0 10]
  cor.scaled <- (abs(data_plot[,which(colnames(data_plot)=="statistic")])-min(abs(data_plot[,which(colnames(data_plot)=="statistic")])))/(max(abs(data_plot[,which(colnames(data_plot)=="statistic")]))-min(abs(data_plot[,which(colnames(data_plot)=="statistic")])))*10
  
  edges <- data.frame(from = data_plot[,which(colnames(data_plot)=="name1")], to = data_plot[,which(colnames(data_plot)=="name2")], value = cor.scaled)
  edges$color <- "black"
  edges$color[edges$value<0] <- "red"
  
  ## plot
  p <- visNetwork(nodes, edges)
  
  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Correlation Network, aes: %s", statname),
      output = list(p)
    )
  ## return
  D
}