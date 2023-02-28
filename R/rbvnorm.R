#' rbvnorm - Gibbs Sampler
#'
#' @param iter number of iterations for gibbs sampler
#' @param mu vector of means
#' @param sigma covariance matrix
#' @importFrom stats rnorm
#'
#' @return a list of the gibbs samples, iterations, mu, and sigma
#' @export
#'
#' @examples rbvnorm(iter = 10000, mu = c(3,6), sigma = matrix(c(8,-2,-2,4), nrow =2, byrow = TRUE))
rbvnorm <- function(iter, mu, sigma){
  x2 <- mu[2]
  df <- matrix(data=NA, nrow=iter,ncol=2)
  colnames(df) <- c("x1", "x2")
  for (i in 1:iter){
    x1 <- rnorm(1, mu[1] + (sigma[1,2]/sigma[2,2])*(x2 - mu[2]), sigma[1,1]-sigma[1,2]*(1/sigma[2,2])*sigma[2,1])
    df[i,1] <- x1
    x2 <- rnorm(1, mu[2] + (sigma[2,1]/sigma[1,1])*(x1 - mu[1]), sigma[2,2]-sigma[2,1]*(1/sigma[1,1])*sigma[1,2])
    df[i,2] <- x2
  }
  invisible(list(gibbs = df, iter = iter, mu = mu, sigma = sigma))
}
