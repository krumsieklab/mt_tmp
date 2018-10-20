# MetaboTools
#
# Activate / deactivate logging
#
# last update: 2018-10-13
# authors: JK
#

require(logging)

# helper function that sends an info message to the "mt" logger
# function also returns the string again, so it can be directly used to store log messages as well
logmsg <- function(msg) { loginfo(msg, logger="mt"); msg }

# main function definition
mt_logging <- function(
  D=NULL,    # SummarizedExperiment, simply being passed through
  console=T  # logging active or inactive?
) {
  
  # set up logging
  logReset()
  if (console) addHandler(writeToConsole, logger = "mt")
  
  # return unchanged object
  if(!is.null(D)) D
  
}



