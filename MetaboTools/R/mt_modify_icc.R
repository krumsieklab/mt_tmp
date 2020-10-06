#' mt_modify_icc
#'
#' Calculate ICC - intraclass correlation
#'
#' @param D \code{SummarizedExperiment} input
#' @param qc_samples Logical expression. Can use fields from \code{colData()}
#' @param col_lab Label for the new icc values to be stored in in rowData
#' @param id_col it is the name of the column in colData containing ids qc_samples
#' @param icc_lmer Optional flag to set lmer flag in ICC function default is FALSE
#' 
#' @return rowData(D)[[col_lab]] contains column with ICC score for each metabolite
#' @examples
#' \dontrun{... %>% mt_modify_icc(qc_samples= (RID%in%RID[duplicated(RID)]), 
#' col_lab = "Dup_icc", id_col='RID') %>% ...}
#'
#' @author RB
#'
#' @importFrom magrittr %>% %<>%
#' @import SummarizedExperiment
#' @importFrom psych ICC
#' @importFrom utils stack
#' @export
mt_modify_icc <- function(D, qc_samples, col_lab, id_col, icc_lmer=F){
  
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(qc_samples)) stop("qc_samples can't be empty!")
  if(missing(id_col)) stop("id_col can't be empty!")
  if(missing(col_lab)) stop("col_lab can't be empty!")
  
  ## APPLY FILTER TO ROW DATA
  qc_samples_q <- dplyr::enquo(qc_samples)
  cd <- colData(D) %>%
    data.frame(row.names = 1:ncol(D)) %>%
    tibble::rownames_to_column("colnames") %>%
    dplyr::filter(!!qc_samples_q)
  
  ## SUBSET SAMPLES
  D1 <- D[, as.numeric(as.matrix(cd$colnames))]
  
  # get the assay of subset with grp column
  qc_data <- bind_cols(D1 %>% assay() %>% t() %>% data.frame(),
                       D1 %>% colData() %>% data.frame() %>% 
                         dplyr::select(!!id_col) %>% data.frame())
  # add ids to the duplicates / multiples
  row_ids <- qc_data %>% 
      dplyr::group_by(!!as.name(id_col)) %>% dplyr::group_rows()
  grp_ids <- data.frame(dup_grp=unlist(lapply(row_ids, FUN=function(x) 1:(length(x))))) %>%
      bind_cols(data.frame(row_id=unlist(row_ids))) %>% .[order(.$row_id), ]
  qc_data$dup_grp <- grp_ids$dup_grp
  # numbers of multiplicates
  num_mults <- qc_data %>% group_by(!!as.name(id_col)) %>% count() 
 # add NAs to the samples if there are differences in number of multiplicates
  if(num_mults %>% pull(n) %>% unique() %>% length() > 1){
    max_mults <- num_mults %>% pull(n) %>% max()
    ids2boot <- num_mults %>% filter(n<max_mults)
    
    addon <- lapply(1:nrow(ids2boot), FUN=function(i){
       tmp <- data.frame(matrix(nrow=(max_mults-ids2boot$n[i]), ncol=ncol(qc_data)))
       names(tmp) <- names(qc_data)
       tmp[['dup_grp']] <- c(1+ids2boot$n[i]):max_mults
       tmp[[id_col]] <- rep(ids2boot[[id_col]][i], nrow(tmp))
       return(tmp)
    }) %>% do.call(rbind, .)
    qc_data <- bind_rows(qc_data, addon)
  }
  # sort by id_col
  qc_data <- qc_data[order(qc_data[[id_col]]), ]
  # names of the metabolites
  mets <- D1 %>% assay() %>% t() %>% data.frame() %>% colnames()
  
  # ICC computation per metabolite
  icc_df <- lapply(mets, function(x){
    icc_met <- qc_data %>% dplyr::select_at(vars('dup_grp', x))
    f <- as.formula(glue::glue(x, ' ~ dup_grp'))
    icc_met <- utils::unstack(icc_met, f)
    icc_val <- tryCatch(psych::ICC(icc_met, missing=F, lmer=icc_lmer)$results[3, 2], silent=T,error = function(err){NA})
    return(icc_val)
  }) %>% do.call(rbind,.) %>% data.frame(icc=.)
  
  # output in rowData
  rowData(D)[[col_lab]] <- icc_df$icc
  
  ## add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Added QC ICC to rowData!")
    )
  ## return
  D
}
