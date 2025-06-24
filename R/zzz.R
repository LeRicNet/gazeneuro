#' @useDynLib gazeneuro, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "gazeneuro package loaded.",
    "\nFor help, see ?gazeneuro or visit https://github.com/lericnet/gazeneuro"
  )
}
