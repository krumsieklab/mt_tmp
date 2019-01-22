# MetaboTools
#
# Activate / deactivate logging
#
# last update: 2018-10-13
# authors: JK
#
# Usage of logger:
#
# Activate it once, either in a SummarizedExperiment pipe, or with a separate call.
# ... %>% mt_logging() %>% ...
# or
# mt_logging()
#
# Then there are three functions for three logging levels:
# mti_logmsg()     - logger "mt", one single message at the end of each mt_ function
# mti_logstatus()  - logger "mts", for internal status updates of the function
# mti_logwarning()    - logger "mtw", for warnings

require(logging)

# helper function that sends an info message to the "mt" logger
# function also returns the string again, so it can be directly used to store log messages as well
mti_logmsg <- function(msg) { loginfo(mti_escape_percent(msg), logger="mt"); msg }
mti_logstatus <- function(msg) { loginfo(mti_escape_percent(msg), logger="mts"); msg }
mti_logwarning <- function(msg) { loginfo(mti_escape_percent(sprintf("WARNING: %s",msg)), logger="mtw"); msg }

# produce sprintf-safe copy of string (create for new loginfo() behaviour, that interprets strings as formats)
mti_escape_percent <- function(txt) gsub('%','%%',txt)

# main function definition
mt_logging <- function(
  D=NULL,    # SummarizedExperiment, simply being passed through
  console=T  # logging active or inactive?
) {
  
  # set up logging
  logReset()
  if (console) addHandler(writeToConsole, logger = "mt")
  if (console) addHandler(writeToConsole, logger = "mts")
  if (console) addHandler(writeToConsole, logger = "mtw")
  
  # return unchanged object
  if(!is.null(D)) D
  
}



