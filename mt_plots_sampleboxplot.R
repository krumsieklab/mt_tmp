# MetaboTools
#
# Boxplot, one box per sample.
# Can be colored by factor.
#
# last update: 2018-10-13
# authors: JK
#

# todo
# - ensure that colorby is a factor or vector of strings

### dependencies
require(reshape2)
require(ggplot2)

source(codes.makepath("packages/metabotools/mt_internal_helpers.R"))

fixorder = function(x){o= unique(as.character(x)); gdata::reorder.factor(x, new.order=o)} # fix order of a factor


### function definition
mt_plots_sampleboxplot <- function(
  D, 
  colorby=NA, 
  title="",
  legend=T,
  ylabel = "Metabolite concentrations"
) {
  
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # sample names with fixed order
  X = t(assay(D))
  samplenames=fixorder(rownames(X))
  
  # prepare for ggplotting
  if (is.na(colorby)) {
    # no coloring
    Xe = cbind(X, data.frame(samplenames=samplenames))
    dfmolten = melt(Xe, id.vars=c("samplenames"))
  } else {
    # with coloring
    Xe = cbind(X, data.frame(groups=get_sample_anno(D,colorby,requireFactor=T), samplenames=samplenames))
    dfmolten = melt(Xe, id.vars=c("groups","samplenames"))
  }
  
  # determine ylim
  stats = sapply(apply(X, 1, boxplot.stats), function(x)x$stats)
  minv = min(stats[1,])
  maxv = max(stats[5,])
  r = maxv-minv
  ylim = c(minv-r/50,maxv+r/50)
  
  # grouped by run
  if (is.na(colorby)) {
    p<-ggplot(data = dfmolten, aes(x=samplenames, y=value)) +
      geom_boxplot(outlier.shape = NA) 
  } else {
    p<-ggplot(data = dfmolten, aes(x=samplenames, y=value)) +
      geom_boxplot(aes(fill=groups),outlier.shape = NA) 
  }
  # visuals
  p <- p + coord_cartesian(ylim=ylim) +
    ggtitle(title) +
    xlab("samples") +
    ylab(ylabel) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.title=element_blank())
  # remove legend?
  if (!legend) p = p + theme(legend.position="none")
  
  # add to metadata plots and return
  metadata(D)$plots %<>% add_to_list(p)
  D
  
}