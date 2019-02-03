# filepath: path of file to be previewed
# filename: name of file to be previewed
preview_roxygen <- function(filepath, filename, rm.temp=T){
  require(roxygen2)
  require(Rdpack)
  
  # create temp folder
  dirpath = tempdir()
  
  # create random string for temp folder
  dirpath = file.path(dirpath, paste0("rp42_",paste0(sample(letters, 10),collapse = "")))
  
  dir.create(dirpath, showWarnings = FALSE)
  dir.create(file.path(dirpath, "man"), showWarnings = FALSE)
  
  # copy the file to be previewed
  file.copy(file.path(filepath,filename),to = dirpath)
  
  mydir = dirpath
  myfiles <- filename
  
  # get parsed source into roxygen-friendly format
  env <- new.env(parent = globalenv())
  rfiles <- sapply(myfiles, function(f) file.path(mydir,f))
  blocks <- unlist(lapply(rfiles, roxygen2:::parse_file, env=env), recursive=FALSE)
  blocks <- lapply(blocks, function(x){
    if(is.null(x$name)) x$name = gsub(pattern = "\\.","_" ,as.character(attributes(x)$call)[2])
    x
  })
  
  # parse roxygen comments into rd files and output then into the "./man" directory
  roc <- roxygen2:::rd_roclet()
  results <- roxygen2:::roclet_process(roc,blocks,env, mydir)
  roxygen2:::roclet_output(roc, results, mydir, options=list(wrap=FALSE), check = FALSE)
  
  options(warn=-1)
  # view html from Rd
  re <-
    sapply(list.files(file.path(dirpath, "man")), function(rd_name)
      Rdpack::viewRd(file.path(dirpath, "man",rd_name),type = "html"))
  options(warn=0)
  
  # remove temporary folder
  if(rm.temp) unlink(dirpath,recursive = T)
  
  if(length(re)==0) return(F)
  
  sapply(re, is.null)
}