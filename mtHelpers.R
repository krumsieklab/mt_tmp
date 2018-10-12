# # MetaboTools
#
# Helper functions.
#
# last update: 2018-10-11
# authors: JK
#

# safely add to end of list, even if list does not exist
addToList <- function(lst, element) {
  # init if needed
  if (is.null(lst) || length(lst)==0) lst=c()
  # add and return
  lst[[length(lst)+1]] = element
  lst
}
