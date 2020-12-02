
#' Faster multivariate normal density using multiple processor cores
#' @param x matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param mean mean vector, default is rep(0, length = ncol(x)).
#' @param sigma covariance matrix
#' @export
fastdmvnorm <- function(x, mean = rep(0, ncol(x)), sigma) {
  cores <- parallel::detectCores(logical = FALSE)
  return(dmvnrm_arma_mc(x, mean, sigma, cores = cores))
}
