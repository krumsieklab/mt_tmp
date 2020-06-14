library(readxl)
library(stringr)
library(SummarizedExperiment)
library(tidyverse)

#' Load Metabolon-format lipidomics data.
#' 
#' Loads lipidomics data from a Metabolon-format Excel file. Needs to be in the original "Client Data Table" format that they deliver.
#'
#' @param file input Excel file
#' @param sheet_list list with the sheet names to read
#' @param copynansheet if given, which sheet to copy the NA pattern from
#'
#' @return Produces an initial SummarizedExperiment, with assay, colData, rowData, and metadata with first entry
#'
#' @examples
#' D <- 
#'   # load data
#'   mt_files_load_metabolon_lipidomics(codes.makepath("Mt/sampledata.xlsx"), sheet_list = c("Species Concentrations","Fatty Acid Concentrations")) %>%
#'   ...
#' 
#' @author EB
#' 
#' @export
mt_files_load_metabolon_lipidomics <- function(
  file,           # Metabolon xls file
  sheet_list     # list of sheet names or numbers to read
) {
  
  # using readxl package:
  raw <- lapply(sheet_list %>% {names(.)=.;.}, function(x){
    read_excel(path=file, sheet=x, col_names = F)
  })
  
  xx <- lapply(sheet_list %>% {names(.)=.;.}, function(x){

    result <- list()
    
    # find metabolite header row and last metabolite row
    imetheader = which(!is.na(raw[[x]][,1]))[3]
    imetlast = max(which(apply(is.na(raw[[x]]),1,sum)<dim(raw[[x]])[2]))
    # find sample header column and last sample row
    isampheader = min(which(!is.na(raw[[x]][5,])))
    isamplast = max(which(apply(is.na(raw[[x]]),2,sum)<dim(raw[[x]])[1]))
    
    # fix overlapping cell
    overl=gsub("\\s+", " ", str_trim(raw[[x]][imetheader,isampheader]))
    overl=strsplit(overl, " ")[[1]]
    overlmet = overl[2]
    overlsamp = overl[1]
    
    # extract metabolite information
    result$metinfo <- read_excel(path=file, sheet=x, col_names = T,
                                      range = cell_limits(ul = c(imetheader+1, 1),
                                                          lr = c(imetlast+1 , isampheader)))
    result$metinfo <- as.data.frame(result$metinfo)
    
    # convert any spaces in the colnames to underscores
    colnames(result$metinfo) <- gsub(" ", "_", colnames(result$metinfo))
    # fix last one
    colnames(result$metinfo)[ncol(result$metinfo)] = overlmet
    
    # extract sample information
    result$sampleinfo = data.frame(t(raw[[x]][5:imetheader,(isampheader+1):isamplast]), stringsAsFactors = F)
    #colnames(result$sampleinfo) = as.list(raw[1:imetheader-1,isampheader])[[1]] # dirty hack, had something to do with the output format of read_excel
    colnames(result$sampleinfo) = c(as.vector(as.matrix(raw[[x]][6:imetheader-1,isampheader])),overlsamp) # dirty hack, had something to do with the output format of read_excel
    rownames(result$sampleinfo) = c()
    result$sampleinfo %<>% mutate_all(parse_guess)
    # convert any spaces in the colnames to underscores
    colnames(result$sampleinfo) <- gsub(" ", "_", colnames(result$sampleinfo))
    
    # extract data
    result$data <- t(raw[[x]][(imetheader+1):imetlast, (isampheader+1):isamplast]) 
    result$data <- as.data.frame(apply(result$data,2, as.numeric))
    
    # as column names, use "Name", if available
    if ("Name" %in% colnames(result$metinfo)) {
      colnames(result$data) = result$metinfo$Name
    } else {
      colnames(result$data) = c()
    }
    
    # as row names, use "CLIENT_IDENTIFIER", if available 
    if ("CLIENT_IDENTIFIER" %in% colnames(result$sampleinfo)) {
      rownames(result$data) = result$sampleinfo$CLIENT_IDENTIFIER
    } else {
      rownames(result$data) = c()
    }
    # add extra column for later merging
    result$data %<>% as.data.frame
    result$data$mergeby <- rownames(result$data)
    
    # set info flags
    result$info$file = file
    result$info$sheet = x
    
    # # copy NanN from another sheet?
    # if (nchar(copynansheet)>0) {
    #   # recursive call
    #   nandf = parseMetabolonFile(file=file, sheet=copynansheet)
    #   # sanity checks
    #   if (!(all.equal(colnames(result$data), colnames(nandf$data)))==T) {
    #     # figure out which are different
    #     l1 = colnames(result$data)
    #     l2 = colnames(nandf$data)
    #     dummy=sapply(1:length(l1), function(i){if(l1[i]!=l2[i]){fprintf('"%s" vs. "%s"\n',l1[i],l2[i])}})
    #     stop('some metabolites are different between data sheet and NaN sheet');
    #   }
    #   if (!all.equal(rownames(result$data), rownames(nandf$data)))
    #     stop('some sample names are different between data sheet and NaN sheet');
    #   # copy over NaNs
    #   result$data[is.na(nandf$data)]=NA
    #   # set info flag
    #   result$info$copynansheet = copynansheet
    #   
    # }

  # add display name
  result$metinfo$name <- result$metinfo$Name
  # # fix variable names  
  # colnames(result$data) <- result$metinfo$Name
  
  result
  
  })
  
  # merge data and annotations from the different data sheets
  result <- list()
  
  # join all sample annotations
  result$sampleinfo <- xx[[sheet_list[1]]]$sampleinfo
  if(length(sheet_list)>1) {
    for (i in 2:length(sheet_list)) {
      result$sampleinfo <- result$sampleinfo %>% 
        full_join(xx[[sheet_list[i]]]$sampleinfo)
    }
  }
  
  # join all metabolite annotations
  result$metinfo <- xx[[sheet_list[1]]]$metinfo
  if(length(sheet_list)>1) {
    for (i in 2:length(sheet_list)) {
      result$metinfo <- result$metinfo %>% 
        rbind(xx[[sheet_list[i]]]$metinfo)
    }
  }
  
  # join all data
  result$data <- xx[[sheet_list[1]]]$data
  if(length(sheet_list)>1) {
    for (i in 2:length(sheet_list)) {
      result$data <- result$data %>% 
        full_join(xx[[sheet_list[i]]]$data, by = "mergeby")
    }
  }
  rownames(result$data) <- result$data$mergeby
  result$data$mergeby <- NULL
  
  # join all info
  result$info$file <- file
  result$info$sheet <- paste(sheet_list, collapse = ", ")
  
  # generate summarized experiment
  D <- SummarizedExperiment(assay    = t(result$data),
                            colData  = result$sampleinfo,
                            rowData  = result$metinfo,
                            metadata = list(sessionInfo=sessionInfo(), parseInfo=result$info))
  
  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("loaded Metabolon lipidomics file: %s, sheets: %s", basename(file), paste(sheet_list, collapse = ", "))
    )
  
  # return
  D
}