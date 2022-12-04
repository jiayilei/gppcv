# source("utilities.R")
#' Title
#'
#' @param X inputs
#' @param y targets
#' @param k covariance function
#' @param sigma2 noise level
#' @param xt test input, 1 dimensional test input data
#'
#' @return
#'        fs mean function
#'        Vfs variance
#'        logp log marginal likelihood
#' @export
#'
#' @examples
gpr <- function(X, y, k, sigma2, xt){
  # compatibility check
  n = dim(X)[1]
  p = dim(X)[2]

  if (length(y) != n){
    stop("dimension of X and Y does not match")
  }
  if (dim(xt)[2] != p){
    stop("dimension of X and Xt does not match")
  }
  if (sigma2 < 0){
    stop("Sigma2 needs to be positive")
  }

  # standardize inputs
  sX <- scale(X)
  sY <- scale(Y)
  X <-matrix(sX, n, p)
  Y <-matrix(sY, n, 1)
  Xt <- (Xt - attr(s, which ="scaled:center"))/ attr(s, which ="scaled:scale")

  # covariance matrix: covariance evaluated at all pairs of training point
  K = matrix(rep(0,n*n), n,n)
  # covariance vector: covariance between test point and the n training points
  ks = matrix(rep(0,n),n)
  for (i in (1:n)){
    for (j in (i:n)){
      K[i,j] = k(X[i,], X[j,])
    }
    ks[i] = k(X[i,], xt)
  }
  tK = t(K)
  diag(tK) <- 0
  K = tK + K

  # cholesky factorization
  L = chol(K + sigma2 * diag(n))

  # find predictive mean
  alpha = solve(t(L)) %*% (solve(L) %*% y)
  fs = t(ks) %*% alpha

  # find predictive variance
  v = solve(L) %*% ks
  Vfs = k(xt, xt) - crossprod(v)

  # calculate log marginal likelihood
  logp = -1/2 * crossprod(y, alpha) - sum(log(diag(L))) - n/2 * log(2*pi)

  return (fs, Vfs, logp)
}


