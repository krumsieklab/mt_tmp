# MetaboTools
#
# Plots correlation network 
#
# last update: 2018-11-08
# authors: EB
#

## dependencies

require(igraph)
require(GGally)
require(network)
require(sna)
source(paste0(codes.makepath("R/networks"),"/writeYED.R", sep=""))

mt_plots_net <- function(
  D,                               # SummarizedExperiment input
  statname,                        # name of the correlation matrix to plot
  corr_filter = p.value < 0.05,    # filter
  export=FALSE,                     # decide if export the network to yEd 
  filename="network.glm"           # name of the yEd file
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
  
  ## apply filter
  if(!missing(corr_filter)){
    mti_logstatus("filter correlations")
    corr_filter_q <- enquo(corr_filter)
    data_plot <- data_plot %>%
      filter(!!corr_filter_q)
  }
  
  data_plot <- data_plot[data_plot$p.adj<0.5, ]
  
  edge_list <- data_plot[,c(which(colnames(data_plot)=="name1"),which(colnames(data_plot)=="name2"),which(colnames(data_plot)=="statistic"))]
  G <- graph.data.frame(edge_list,directed=FALSE)
  A <- abs(as_adjacency_matrix(G, type="both", names=TRUE, sparse=FALSE, attr="statistic"))
  
  net <- network(A, directed = FALSE, names.eval = "weights")
  network.vertex.names(net) = colnames(A)
  set.edge.attribute(net, "color", ifelse(net %e% "weights" > 0, "black", "red"))
  
  ## create plot
  p <- ggnet2(net, label = TRUE, node.size = 6) #, edge.color = "color" , edge.size = "weights"
  
  ## export network to file
  if(export) {
    sprintf("network exported to file %s", filename)
    A.adj <- A
    A.adj[A.adj!=0] <- 1
    pos <- c(0, 0, 0)
    neg <- c(1, 0, 0)
    col <- rbind(pos, neg)
    writeYED(A = A.adj, outfile = filename,    # not working for weighted adjacency right now, only plots nodes and labels, no edges
             options = set_options(nodelabels=colnames(A)))
  }
  
  
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