library(GGally)
library(sna)
library(ggnetwork)
library(visNetwork)
library(dils)

#' mt_plots_net
#' 
#' Creates a network object
#'
#' @param D \code{SummarizedExperiment} input
#' @param statname name of the test to take correlations from
#' @param corr_filter filter for correlation values to plot
#' @param node_coloring name of the test to use for node coloring
#' @param save.html filename of visnetwork html file. If empty, no html saved
#' 
#' @return assay: not altered
#' @return $result: network ggplot + visnetwork plot
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
  node_coloring,                   # name of the statistical test to use for node coloring
  save.html                        # filename of visnetwork html
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
  
  # define node attributes
  nodes <- as.data.frame(unique(rbind(cbind(ids=data_plot$var1,label=data_plot$name1),cbind(ids=data_plot$var2,label=data_plot$name2))))
  
  if(!(missing(node_coloring))) {
    test <- mti_get_stat_by_name(D, node_coloring);
    test <- test[match(nodes$ids, test$var),];
    nodes$map <- sign(test$statistic)*log10(test$p.value);
    nodes$node_color <- colorRampPalette(c('blue', 'white', 'red'))(length(nodes$map))[rank(nodes$map)]
  } else nodes$node_color <- rep("darkblue",times=length(nodes$ids))

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
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$color <- "black"
  edges$color[edges$value<0] <- "red"
  
  ## plot
  e <- edges
  n <- data.frame(id=nodes$label, label= nodes$label, color=nodes$node_color)
  p_vis <- visNetwork(n,e)
  
  # greate ggnetwork object
  df <- list()
  df$edges <- data.frame(from = data_plot[,which(colnames(data_plot)=="var1")], to = data_plot[,which(colnames(data_plot)=="var2")], weight = cor.scaled/10)
  df$vertices <- nodes$label
  
  adj <- as.matrix(AdjacencyFromEdgelist(df$edges, check.full = TRUE))
  colnames(adj[[1]])<- as.character(nodes$label[match(adj[[2]],nodes$ids)])
  rownames(adj[[1]])<- as.character(nodes$label[match(adj[[2]],nodes$ids)])
  mm.net <- network(adj[[1]], layout = "kamadakawai", directed = FALSE)
  
  test <- mti_get_stat_by_name(D, node_coloring);
  test <- test[match(adj[[2]], test$var),];
  map <- sign(test$statistic)*log10(test$p.value);
  
  mm.col <- c("positive" = "#000000", "negative" = "#0000FF")
  x <- data_plot$statistic
  x[x>=0] <- "positive"
  x[x<0] <- "negative"
  mm.net %e% "pcor" <- abs(data_plot[,which(colnames(data_plot)=="statistic")])
  mm.net %e% "pos" <- x
  mm.net %v% "strength" <- map
  
  # add edge color (positive/negative)
  # add black circle around nodes
  p <- ggplot(mm.net, aes(x, y, xend = xend, yend = yend)) +
    geom_edges(color="black",aes(size=pcor)) +
    geom_nodes(aes(color = strength), size = 7) +
    scale_color_gradient2(low = "#0000FF", mid="white",high = "#FF0000") +
    geom_nodetext(color="grey50",aes(label = vertex.names),
                  size = 3, vjust = -0.6) +
    theme_blank() +
    theme(legend.position = "bottom")

  # if save.html given, save visnetwork to html
  if (!missing(save.html)){
    visSave(graph = p_vis, file = save.html)
  }

  ## add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Correlation Network, aes: %s", statname),
      output = list(p),
      output2 = p_vis
    )
  ## return
  D
}
