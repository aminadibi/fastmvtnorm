binorm_pdf <- function(x, sigma) {
  s1 <- sqrt(sigma[1, 1])
  s2 <- sqrt(sigma[2, 2])
  rho <- sigma[1, 2] / (s1 * s2)
  Z <- x[1] ^ 2 / sigma[1, 1] - 2 * rho * x[1] * x[2] / (s1 * s2) + x[2] ^ 2 / sigma[2, 2]
  densRes <- 1 / (2 * pi * s1 * s2 * sqrt(1 - rho ^ 2)) * exp(- Z / (2 * (1 - rho ^ 2)))
  return(densRes)
}


#' Faster multivariate normal density using multiple processor cores
#' @param x matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param mean mean vector, default is rep(0, length = ncol(x)).
#' @param sigma covariance matrix
#' @export
fastdmvnorm <- function(x, mean = rep(0, length(x)), sigma) {
  if (length(x) == 2) {
    result <- binorm_pdf(x, sigma)
  } else {
    cores <- detectCores(logical = FALSE)
    result <- dmvnrm_arma_mc(x, mean, sigma, cores = cores)
  }

  return(result)
}
