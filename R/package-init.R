#' @useDynLib fastmvtnorm, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("fastmvtnorm", libpath)
}
