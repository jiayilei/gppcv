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
  standardized_out <- standardize(X, y, Xt, yt)

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


# standardize input
standardize <- function(X, y, Xt, yt){
  n = dim(X)[1]
  p = dim(X)[2]

  sX <- scale(X)
  sY <- scale(y)
  X <-matrix(sX, n, p)
  y <-matrix(sY, n, 1)
  Xt <- (Xt - attr(sX, which ="scaled:center"))/ attr(sX, which ="scaled:scale")
  yt <- (Yt - attr(sY, which ="scaled:center"))/ attr(sY, which ="scaled:scale")

  return (list(X = X, y = y, Xt = Xt, yt = yt))
}


# returns Euclidean distance between two input data
get_r <- function(x1, x2){
  return (sqrt(sum((x1- x2)^2)))
}


# Square Exponential Kernel
se_kernel <- function(l, r = 1){
  fun <- function(r){
    return (exp(-0.5 * (r/l)^2))
  }
  return (fun)
}


# Matern kernel
matern_kernel <- function(l, v, r = 1){
  fun <- function (r){
    left <-  1 / gamma(v) / 2^(v-1)
    mid <- (sqrt(2*v)/ l * r)^ v
    right <- besselK(sqrt(2*v) * r / l, nu = v)
    return (left * mid * right)
  }
  return (fun)
}


# exponential kernel
exp_kernel <- function(l, r=1){
  fun <- function(r){
    return (exp(-r/l))
  }
  return (fun)
}


# pick the kernel
pick_kernel <- function (method = c('se', 'm', 'exp'), para_list){
  # if (method == 'se'){
  #   return (se_kernel)
  # }
  # if (method == 'm'){
  #   return (matern_kernel)
  # }
  # if (method == 'exp'){
  #   return (exp_kernel)
  # }

  if (method == 'se'){
    return (do.call(se_kernel, para_list))
  }
  if (method == 'm'){
    return (do.call(matern_kernel, para_list))
  }
  if (method == 'exp'){
    return (do.call(exp_kernel, para_list))
  }

}
