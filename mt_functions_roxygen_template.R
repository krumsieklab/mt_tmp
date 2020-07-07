# Roxygen template for functions - KEEP THIS FILE OUTSIDE PACKAGE

# exported functions
# Notes on the exported function template below:
#   1. import and importFrom are included as needed
#   2. @export is required for function to be accessible to package users
#   3. put example code in between braces of \dontrun{} if the code cannot be run on its own - e.g. if following ... or requires data not automatically loaded
# Exported functions template below the line
# ------------------------------
#' Title / Header of Function
#'
#' Description of function
#'
#' @param function_argument Description of function argument. Add one of these for each argument
#'
#' @return what the function returns
#'
#' @examples
#' # example of how to run function
#' \dontrun{}
#'
#' @author author_initials
#'
#' @import name_of_package
#' @importFrom name_of_package function_or_operator
#'
#' @export


# internal functions
# Notes on the internal function template below:
#   1. import and importFrom are included as needed
#   2. @noRd indicates that documentation will not be written for this function
#   3. @examples not included - not necessary as won't be accessible to user
# internal functions template below the line
# ---------------------------
#' Title / Header of Function
#'
#' Description of function
#'
#' @param function_argument Description of function argument. Add one of these for each argument
#'
#' @return what the function returns
#'
#' @author author_initials
#'
#' @import name_of_package
#' @importFrom name_of_package function_or_operator
#'
#' @noRd
