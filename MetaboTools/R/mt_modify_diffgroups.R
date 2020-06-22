#' Generate a new data matrix that is one group minus a second group. (e.g. difference between two timepoints).
#'
#' Takes the name of a sample factor variable (from colData), and two factor levels. Creates a new data matrix that
#' is made of the second minus the first levels. Will only produce results for cases where both
#'
#' @param D \code{SummarizedExperiment} input
#' @param idvar sample ID variable  (e.g. patient ID)
#' @param groupvar sample variable from colData (e.g. timepoint)
#' @param grp1 group level 1 (e.g. "day1")
#' @param grp2 group level 2 (e.g. "day2")
#'
#' @export

mt_modify_diffgroups <- function(
  D,
  idvar,
  groupvar,
  grp1,
  grp2
) {

  # validate arguments
  stopifnot("SummarizedExperiment" %in% class(D))

  # if data is logged, then subtract, otherwise divide
  if (length(mti_res_get_path(D, c("pre","trans","log"))) > 0){
    op <- `-`
    opstr <- "-"
  } else {
    op <- `/`
    opstr <- "/"
  }


  # shortcuts
  cd = colData(D)
  X = assay(D)

  # verify that all fields exist
  if (!(idvar %in% colnames(cd))) stop(sprintf("'%s' not found in sample annotations.", idvar))
  if (!(groupvar %in% colnames(cd))) stop(sprintf("'%s' not found in sample annotations.", groupvar))

  # find matching samples
  s1 = cd[[groupvar]] == grp1
  s2 = cd[[groupvar]] == grp2

  # extract first group
  ids1 = cd[[idvar]][s1]
  dat1 = X[,s1]
  # extract second group
  ids2 = cd[[idvar]][s2]
  dat2 = X[,s2]

  # map the two groups
  idc = intersect(ids1,ids2)
  m1 = match(idc,ids1)
  m2 = match(idc,ids2)

  # subtract or divide data frames
  diffmat = op(dat2[,m2], dat1[,m1])
  # copy the clinical data of the first group 1
  cd_sub <- cd[s1,]
  cd_new <- cd_sub[match(idc, cd_sub[[idvar]]),]

  # create new SummarizedExperiment object
  D <- SummarizedExperiment(assay    = diffmat,
                            rowData  = rowData(D),
                            colData  = cd_new,
                            metadata = metadata(D))

  # add status information
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = glue('created new dataset, id var: {idvar}, subtracting {groupvar}, group {grp2}{opstr}{grp1}')
    )

  # return
  D


}
