# MetaboTools
#
# Produce equalizer plot.
# Requires two nested annotations, e.g. SuperPathway and SubPathway, or SubPathway and Metabolites
#
# last update: 2019-01-03
# authors: JK
#

# TODO: legend labels
# JAN TODO NEXT: MAKE WORK WITH SUB AND METABOLITES, WORKSPACE SAVED

library(purrr)
library(glue)

mt_plots_equalizer <- function(
  D1,       # SummarizedExperiment input 1
  comp1,    # name of first comparison output to take arguments from     first one has to be the less granular one (e.g. D1 super, D2 sub)
  D2,       # SummarizedExperiment input 2
  comp2,    # name of second comparison output to take arguments from
  th = 2,   # log10(p.value) threholds for red dashed lines
  clrs = c("#9494FF","red") # colors for sub and super pathways
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D1))
  stopifnot("SummarizedExperiment" %in% class(D2))
  stopifnot(comp1!=comp2)
    
  # get results
  res1 <- mti_get_stat_by_name(D1, comp1) 
  res2 <- mti_get_stat_by_name(D2, comp2) 
  
  # shortcuts
  rd1 = rowData(D1)
  rd2 = rowData(D2)
  
  # find field in first that contains all the second (e.g. all subpathways in the metabolite rowData)
  col2 <- rd2 %>% sapply(function(x) all((x %>% na.omit()) %in% rd1[["name"]]) )
  if (sum(col2)>1) stop(sprintf("Multiple columns in first rowData map to names of second: %s", paste0(colnames(rd2)[col2],collapse=", ")))
  if (sum(col2)==0) stop(sprintf("No columns in first rowData map to names of second. Specified correct rowData frames?"))
  
  
  ##### Mustafa's code starts here
  
  # sub pathways stats
  df = data.frame(res2[,c("statistic", "p.value")], rd2[,-1])
  # super pathways stats
  df2 = data.frame(SUPER_PATHWAY=gsub("\\.", " ", res1$var),res1[,-(1:3)], stringsAsFactors = F)
  
  # create x-axis
  df$x = abs(log10(df$p.value)) * sign(df$statistic)
  # x axis limits
  a = max(abs(df$x))
  
  # main facetted plot
  gg<-
    ggplot(df, aes(x = x, y = SUB_PATHWAY)) +
    geom_vline(xintercept = 0, color ="gray") +
    geom_vline(xintercept = c(-th,th), color ="tomato", lty = 2) +
    geom_point(pch = 22, fill = clrs[1], size = 3) +
    facet_grid(SUPER_PATHWAY~. , scales = "free_y", space = "free_y") +
    theme(strip.background =element_rect(fill=NA), 
          strip.text = element_text(colour = 'black', face = "bold"), 
          strip.text.y = element_text(angle = 0, hjust = 0), 
          panel.grid.major.y = element_line(color ="gray"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_rect(fill=NA, color ="black")) + 
    scale_x_continuous(limits = c(-a,a))
  
  # add super pathways 
  df.super = dplyr::summarise(dplyr::group_by(df,SUPER_PATHWAY), 
                              yy = mean(as.numeric(factor(SUB_PATHWAY))),
                              xx = mean(x))
  if(!is.null(df2)){
    df.super = dplyr::inner_join(df.super, df2, "SUPER_PATHWAY")
    # create x-axis
    df.super$xx = abs(log10(df.super$p.value)) * sign(df.super$statistic)
  }
  
  gg = gg + geom_point(data = df.super, aes(x = xx,y = yy),pch = 22, 
                       fill = clrs[2], size = 5, alpha = 0.7) 
  
  # add legend
  df.legend = data.frame(x = rep(1,2), y = rep(NA,2), class = c("Super pathways", "Sub pathways"))
  p <- gg + geom_point(data = df.legend,aes(x= x,y=y, fill =class),pch = 22, size = 4) + 
    labs(fill="", x = expression(paste("directed |log"[10],"p|")), y = "") + 
    scale_fill_manual(values = clrs) + 
    theme(legend.position = "top", legend.key = element_blank(), 
          legend.direction = "vertical", legend.justification = c(0,0))
  
  
  # add status information & plot
  funargs <- mti_funargs()
  metadata(D1)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("equalizer plot of '%s' and '%s'", comp1, comp2),
      output = list(p)
    )
  
  # return
  D1
 
  
}









