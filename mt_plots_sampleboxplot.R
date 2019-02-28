# MetaboTools
#
# Boxplot, one box per sample.
# Can be colored by factor.
#
# last update: 2018-10-13
# authors: JK, JZ
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
  D,         # SummarizedExperiment input
  plottitle="Sample boxplot",  # title of boxplot
  legend=T,  # keep legend?  [could be removed]
  ylabel = "Metabolite concentrations",  # y axis label
  logged=F,  # plot logged
  ...        # additional arguments directly passed to aes() of ggplot
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # logged?
  Dplot = D
  if (logged) {
    assay(Dplot) <- log2(assay(Dplot))
    ylabel = sprintf("%s [log2]", ylabel)
  }
  
  # generate ggplot
  p <- Dplot %>%
    mti_format_se_samplewise() %>%
    gather(metab, value, one_of(rownames(D))) %>%
    ggplot(aes(x = primary, y = value, ...)) +
    geom_boxplot() +
    ylab(ylabel) +
    ggtitle(plottitle) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # todo add rowname 
  
  # remove legend?
  if (!legend) p = p + theme(legend.position="none")
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("sample boxplot, aes: %s", mti_dots_to_str(...)),
      output = list(p)
    )
  
  # return
  D
  
  
}
