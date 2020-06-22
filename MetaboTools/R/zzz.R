.onLoad <- function(libname, pkgname){
  if(exists("mt_files_data_xls", envir = globalenv())){  # reference function to check for previously loaded mt version
    stop("ANOTHER VERSION OF METABOTOOLS WAS DETECTED IN THE GLOBAL ENVIRONMENT.\nPlease completely detach any other versions
of MetaboTools from the global environment before reloading library.")
  }
}
