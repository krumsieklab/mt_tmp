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
  title="",  # title of boxplot
  legend=T,  # keep legend?  [could be removed]
  ylabel = "Metabolite concentrations",  # y axis label
  ...        # additional arguments directly passed to aes() of ggplot
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  # generate ggplot
  p <- D %>%
    mti_format_se_samplewise() %>%
    gather(metab, value, one_of(rownames(D))) %>%
    ggplot(aes(x = primary, y = value, ...)) +
    geom_boxplot()
  # todo add rowname 
  
  # remove legend?
  if (!legend) p = p + theme(legend.position="none")
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("sample boxplot, aes: %s", mti_dots_to_str(...)),
      output = p
    )
  
  # return
  D
  
  
}
