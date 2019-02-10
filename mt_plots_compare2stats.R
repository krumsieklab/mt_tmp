library(ggrepel)

mt_plots_compare2stats <- function(
  D1,
  stat1,
  filter1,
  stat2,
  filter2,
  filterop="AND",
  D2 = D1,
  return.plot.only=F  # return only the plot object. note: setting this to true makes the function non-MT pipeline compatible.
) {
  
  ## check input
  stopifnot("SummarizedExperiment" %in% class(D1))
  stopifnot("SummarizedExperiment" %in% class(D2))
  filter1q <- enquo(filter1)
  filter2q <- enquo(filter2)
  
  if (!(filterop %in% c("AND","OR"))) stop("filterop must be 'AND' or 'OR'")
  
  ## obtain the two stats structures
  s1 <- mti_get_stat_by_name(D1, stat1, fullstruct=T)
  s1t <- s1$table
  s2 <- mti_get_stat_by_name(D2, stat2, fullstruct=T)
  s2t <- s2$table
  
  ## add directed p-value, filter, and merge
  s1t$dp1 <- -log10(s1t$p.value) * sign(s1t$statistic)
  s1t$filtered1 <- s1t$var %in% (s1t %>% filter(!!filter1q))$var
  s2t$dp2 <- -log10(s2t$p.value) * sign(s2t$statistic)
  s2t$filtered2 <- s2t$var %in% (s2t %>% filter(!!filter2q))$var
  st <- merge(s1t, s2t, by='var')
  
  
  st <- merge(st, rowData(D1), by.x="var", by.y='row.names', all.x=T) # add names
  
  
  # combine filters
  if (filterop=="AND")
    # AND
    st$filtered = as.numeric(st$filtered1 & st$filtered2)
  else if (filterop=="OR")
    # OR
    st$filtered = as.numeric(st$filtered1) + as.numeric(st$filtered2)
  else
    stop("bug")
  
  ## create axis labels
  if ("groups" %in% names(s1) && length(s1$groups)==2) {
    xlabel = sprintf("%s high <--   dir. log10(p-value)   --> %s high", s1$groups[1], s1$groups[2])
  } else {
    xlabel = 'directed log10(p-value)'
  }
  if ("groups" %in% names(s2) && length(s2$groups)==2) {
    ylabel = sprintf("%s high <--   dir. log10(p-value)   --> %s high", s2$groups[1], s2$groups[2])
  } else {
    ylabel = 'directed log10(p-value)'
  }
  
  
  ## plot
  st <- as.data.frame(st)
  p <- st %>% 
    ggplot(aes(x=dp1,y=dp2,color=as.factor(filtered))) + 
    geom_point() + 
    labs(color='filtered') +
    geom_text_repel(data=filter(st, filtered>0), aes(label=name), size=3, colour = "black") + 
    xlab(xlabel) + ylab(ylabel)
  
  
  if (!return.plot.only) {
    ## add status information & plot
    funargs <- mti_funargs()
    metadata(D1)$results %<>% 
      mti_generate_result(
        funargs = funargs,
        logtxt = sprintf("comparison plot between '%s' and '%s'", stat1, stat2),
        output = list(p),
        output2 = st
      )
    ## return
    D1
  } else {
    p
  }
  
}
