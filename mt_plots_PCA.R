# MetaboTools
#
# PCA plot, show first two PCs plotted against each other.
# - Coloring by factor or numerical.
# - Symbols by factor.
#
# last update: 2018-10-13
# authors: JK
#


require(ggplot2)

mt_plots_PCA <- function(
  D,
  colorby=NA,
  shapeby=NA,
  title=""
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # extract data and verify that there are no NA values
  X = t(assay(D))
  if (any(is.na(X))) stop("Data matrix for PCA cannot contain NAs")
  
  # PCA
  pca <- prcomp(x=as.matrix(X)) 
  expvar = (pca$sdev)^2 / sum(pca$sdev^2) 
  # assemble data frame
  df = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
  
  # prep
  pc1name = sprintf('PC1 (%.1f%%)', expvar[1]*100)
  pc2name = sprintf('PC2 (%.1f%%)', expvar[2]*100)
  
  # build aesthetics
  a <- aes(x=PC1,y=PC2)
  # add coloring?
  if (!is.na(colorby)) {
    a <- c(a, aes(color=color))
    class(a) <- "uneval"
    df$color <- get_sample_anno(D,colorby,requireFactor=F)
  } 
  # add shapeby?
  if (!is.na(shapeby)) {
    a <- c(a, aes(shape=shape))
    class(a) <- "uneval"
    df$shape <- get_sample_anno(D,shapeby,requireFactor=F)
  } 
  
  # class(a.merged) <- "uneval"
  
  # plot
  p <- ggplot(data=df)
  p = p + geom_point(size=3,mapping=a)
  p = p + xlab(pc1name) + ylab(pc2name) + ggtitle(title) + theme(legend.title=element_blank())
  
  # add to metadata plots and return
  metadata(D)$plots %<>% add_to_list(p)
  D
  
  # if (!useshapes) {
  #   df = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], grps=grps)
  #   p<-ggplot(data=df) + 
  #     geom_point(size=4,aes(x=PC1,y=PC2,fill=grps), pch=21) + 
  #     xlab(pc1name) + ylab(pc2name) + ggtitle(title)
  # } else {
  #   df = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], grps=factor(grps), grps2 = as.factor(grps2))
  #   p<-ggplot(data=df, aes(x=PC1,y=PC2,shape=grps2)) + 
  #     geom_point(size=4,aes(colour=grps)) + 
  #     xlab(pc1name) + ylab(pc2name) + ggtitle(title)
  # }
  # # pairing lines?
  # if (!(is.na(pairing)) && length(pairing)>0) {
  #   dflines <- data.frame(x=pca$x[,1], y=pca$x[,2],grps=pairing)
  #   p = p+geom_line( data=dflines, aes(x=x,y=y,group=grps), alpha=0.7, colour="#999999")
  # }
  
  
  

  
  
}
