# Logging control codes. If not explicitely called, MT will perform standard logging to the console.
library(logging)

# helper function that sends an info message to the "mt" logger
# function also returns the string again, so it can be directly used to store log messages as well
mtm_logmsg <- function(msg) { loginfo(mti_escape_percent(msg), logger="mt"); msg }
mtm_logstatus <- function(msg) { loginfo(mti_escape_percent(msg), logger="mts"); msg }
mtm_logwarning <- function(msg) { loginfo(mti_escape_percent(sprintf("WARNING: %s",msg)), logger="mtw"); msg }

# produce sprintf-safe copy of string (create for new loginfo() behaviour, that interprets strings as formats)
mti_escape_percent <- function(txt) gsub('%','%%',txt)

# main function definition
mtm_logging <- function(
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



