library(purrr)
library(glue)

#' 'Equalizer' plots.
#' 
#' Creates a nested plot based on metabolite annotations, e.g. of super-/ and sub-pathways, or of sub-pathways and metabolites.
#'
#' @param D1 \code{SummarizedExperiment} input 1, the coarse one
#' @param comp1  name of first comparison output to take arguments from, the coarse one [first one has to be the less granular one (e.g. D1 super, D2 sub)]
#' @param D2 \code{SummarizedExperiment} input 2, the fine one
#' @param comp2 name of second comparison output to take arguments from, the fine one
#' @param legend.fine fine label to be plotted
#' @param legend.coarse coarse legend to be plotted
#' @param vertline.fine filter expression where to draw the red, dashed line, for fine. default: 0.05
#' @param vertline.coarse filter expression where to draw the red, dashed line, for coarse. default: 0.05
#' @param clrs colors for fine and coarse, default: c("#9494FF","red") (light blue and red)
#'
#' @return $result: plot, equalizer
#' 
#' @examples
#' # super-pathway / sub-pathway equalizer
#' # sub-pathway analysis must already be stored in D_sub, and this is part of the super-pathway pipeline, with a result already in 'comp'
#'  ... %>%
#'  mt_plots_equalizer(
#'   comp1='comp', 
#'   D2=D_sub, 
#'   comp2=='comp',
#'   legend.fine="sub pathway", 
#'   legend.coarse='super pathway',
#'   vertline.fine = p.adj < 0.1,
#'   vertline.coarse = p.adj < 0.1) %>% 
#' ...
#' 
#' @author JK, MB
#' 
#' @export
mt_plots_equalizer <- function(
  D1,       # SummarizedExperiment input 1, the coarse one
  comp1,    # name of first comparison output to take arguments from, the coarse one [first one has to be the less granular one (e.g. D1 super, D2 sub)]
  D2,       # SummarizedExperiment input 2, the fine one
  comp2,    # name of second comparison output to take arguments from, the fine one
  legend.fine, # fine label to be plotted
  legend.coarse = NULL, # coarse legend to be plotted
  vertline.fine = p.adj < 0.05, # filter expression where to draw the red, dashed line, for fine
  vertline.coarse = p.adj < 0.05, # filter expression where to draw the red, dashed line, for coarse
  clrs = c("#9494FF","red") # colors for sub and super pathways
) {
  
  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D1))
  stopifnot("SummarizedExperiment" %in% class(D2))
  # stopifnot(comp1!=comp2)
  
  # get results
  res1 <- mti_get_stat_by_name(D1, comp1) 
  res2 <- mti_get_stat_by_name(D2, comp2) 
  
  ##### Mustafa's code starts here
  
  # rm "name" duplicates
  frm_ <- function(x){
    xnames = setdiff(colnames(x),"name")
    x[,c("name", xnames[!sapply(x[,xnames,drop =F],identical, x$name)]),drop =F]
  }
  
  # shortcuts
  rd1 = rowData(D1) %>% as.data.frame %>% frm_
  rd2 = rowData(D2) %>% as.data.frame %>% frm_
  
  
  # find field in first that contains all the second (e.g. all subpathways in the metabolite rowData)
  col2 <- rd2 %>% sapply(function(x) all((x %>% na.omit()) %in% rd1[["name"]]) ) 
  
  if (sum(col2)>1) stop(sprintf("Multiple columns in first rowData map to names of second: %s", paste0(colnames(rd2)[col2],collapse=", ")))
  if (sum(col2)==0) stop(sprintf("No columns in first rowData map to names of second. Specified correct rowData frames?"))
  
  
  col2 = col2 %>% which %>% names
  colnames(rd2)[colnames(rd2) == col2] = "COARSE"
  # if legend.coarse not given 
  if(is.null(legend.coarse)) legend.coarse = col2
  
  colnames(rd1)[1] = "COARSE"
  colnames(rd2)[1] = "FINE" 
  
 
  
  # df: data frame includes columns: "SUB_PATHWAY", "SUPER_PATHWAY", "statistic", "p.value"
  # name.df: primary key(column) name in df
  # df2: data.frame for super pathways, columns: SUPER_PATHWAY", "statistic", "p.value"
  # name.df2: foreign key(column) name in df2
  # th: log10(p.value) threholds for red dashed lines
  # clrs: colors for sub and super pathways
  
  mti_plot_equalizer_gg <- function(df, name.df="SUB", df2=NULL, name.df2="SUPER" ){
    
    # create x-axis
    df$x = abs(log10(df$p.value)) * sign(df$statistic)
    # x axis limits
    a = max(abs(df$x))
    
    # find x coordinates for cutoff lines
    xfine <- res2 %>% filter(!!enquo(vertline.fine)) %>% .$p.value %>% max()
    xcoarse <- res1 %>% filter(!!enquo(vertline.coarse)) %>% .$p.value %>% max()
    
    # main facetted plot
    gg<-
      ggplot(df, aes(x = x, y = FINE)) +
      geom_vline(xintercept = 0, color ="gray") +
      geom_vline(xintercept = c(-log10(xfine),log10(xfine)), color=clrs[1], alpha=0.4) + 
      geom_vline(xintercept = c(-log10(xcoarse),log10(xcoarse)), color=clrs[2], alpha=0.4) +
      geom_point(pch = 22, fill = clrs[1], size = 3) +
      facet_grid(COARSE~. , scales = "free_y", space = "free_y") +
      theme(strip.background =element_rect(fill=NA), 
            strip.text = element_text(colour = 'black', face = "bold"), 
            strip.text.y = element_text(angle = 0, hjust = 0), 
            panel.grid.major.y = element_line(color ="gray"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.background = element_rect(fill=NA, color ="black")) + 
      scale_x_continuous(limits = c(-a,a))
    
    # add super pathways 
    df.super = dplyr::summarise(dplyr::group_by(df,COARSE), 
                                yy = mean(as.numeric(factor(FINE))),
                                xx = mean(x))
    
    if(!is.null(df2)){
      # if NA exist in super but in df2
      if(any(is.na(df.super$COARSE)) & !any(is.na(df2$COARSE))){
        df2 = rbind(df2, rep(NA,ncol(df2)))
      }
      
      df.super = dplyr::inner_join(df.super, df2, "COARSE")
      # create x-axis
      df.super$xx = abs(log10(df.super$p.value)) * sign(df.super$statistic)
    }
    
    
    gg = gg + geom_point(data = df.super, aes(x = xx,y = yy),pch = 22, 
                         fill = clrs[2], size = 5, alpha = 0.7) 
    
    # add legend
    df.legend = data.frame(x = rep(1,2), y = rep(NA,2), class = factor(c(name.df2, name.df),levels = c(name.df, name.df2)))
    gg + geom_point(data = df.legend,aes(x= x,y=y, fill =class),pch = 22, size = 4) + 
      labs(fill="", x = expression(paste("directed log10(p)")), y = "") + 
      scale_fill_manual(values = clrs) + 
      theme(legend.position = "top", legend.key = element_blank(), 
            legend.direction = "vertical", legend.justification = c(0,0))
    
  }
  
  
  
  p =  mti_plot_equalizer_gg(df = data.frame(rd2, res2), name.df = legend.fine,
                             df2 = data.frame(rd1, res1), name.df2 = legend.coarse )
  
  ## ADD AXIS GROUPS
  d <- mti_get_stat_by_name(D1, comp1, fullstruct=T)
  if ("groups" %in% names(d) && length(d$groups)==2) {
    xlabel <- sprintf("%s high <--     directed log10(p)     --> %s high", d$groups[1], d$groups[2])
    p <- p + xlab(xlabel)
  }
  

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









