# source("utilities.R")
#' Title
#'
#' @param X inputs
#' @param y targets
#' @param k covariance function
#' @param sigma2 noise level
#' @param xt test inputs
#' @param yt test targets
#'
#' @return
#'        fs mean function
#'        Vfs variance
#'        logp log marginal likelihood
#' @export
#'
#' @examples
gpr_standardized <- function(X, y, k, sigma2, xt, yt, K, ks){
  # compatibility check
  n = dim(X)[1]
  p = dim(X)[2]
  nt = dim(Xt)[1]

  # cholesky factorization
  L = chol(K + sigma2 * diag(n))

  # find predictive mean
  alpha = solve(t(L)) %*% (solve(L) %*% y)
  fs = ks %*% alpha

  # find predictive variance
  v = solve(L) %*% t(ks)
  r = matrix(rep(0, nt), nt)
  for (i in 1:nt){
    r[i] = get_r(xt[i,], xt[i,])
  }
  Vfs = k(r) - colSums(v^2)

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
  yt <- (yt - attr(sY, which ="scaled:center"))/ attr(sY, which ="scaled:scale")

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
  # todo: check dimension of inputs

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

# get covariance_matrix of K and Ks
covariance_mats <- function(X, Xt, k){
  n = dim(X)[1]
  p = dim(X)[2]
  nt = dim(Xt)[1]

  # covariance matrix: covariance evaluated at all pairs of training point
  K = matrix(rep(0,n*n), n,n)
  # covariance vector: covariance between test point and the n training points
  ks = matrix(rep(0,nt * n), nt, n)
  for (i in (1:n)){
    for (j in (1:n)){
      r = get_r(X[i,], X[j,])
      K[i,j] = k(r)
    }
    for (k in (1: nt)){
      r = get_r(X[i,], xt[k,])
      ks[k,i] = k(r)
    }

  }
  tK = t(K)
  diag(tK) <- 0
  K = tK + K

  return (list(K = K, ks = ks))
}


# sequence of k's to fit the model and compare different kernels performance
gpr_seq_kernels <-function(X, y, k, sigma2, Xt, yt){
  # compatibility check
  n = dim(X)[1]
  p = dim(X)[2]
  nt = dim(Xt)[1]

  if (length(y) != n){
    stop("dimension of X and Y does not match")
  }
  if (dim(Xt)[2] != p){
    stop("dimension of X and Xt does not match")
  }
  if (sigma2 < 0){
    stop("Sigma2 needs to be positive")
  }

  # standardize inputs
  standardized_out <- standardize(X, y, Xt, yt)
  X = standardized_out$X
  y = standardized_out$y
  Xt = standardized_out$Xt
  yt = standardized_out$yt

  # get K and Ks, covariance matrix
  # browser()


  # standardized gaussian process regression for each k
  mse <- matrix(rep(0, length(k)))
  for (i in 1 : length(k)){
    cov_out < covariance_mats(X, Xt, k[i])
    K = cov_out$K
    ks = cov_out$ks
    gpr_out <- gpr_standardized(X, y, k[[i]], sigma2, Xt, yt, K, ks)
    # evaluate mse for each kernels
    mse[i] = sum((gpr_out$fs - yt)^2)/nt

  }
  return (list(k=k, mse = mse))

}

# cv
gpr_cv <- function(X, y, k, sigma2, num_folds){
  n = dim(X)[1]
  fold_ids <- sample((1:n) %% num_folds + 1, n)

  nkernels <- length(k)
  cv_folds <- matrix(rep(num_folds * nkernels), num_folds, nkernels)
  for (fold in 1:num_folds) {
    Xtrain <- X[fold_ids != fold, ]
    ytrain <- y[fold_ids != fold]

    Xtest <- X[fold_ids == fold, ]
    ytest <- y[fold_ids == fold]

    # browser()
    gpr_seq_out <- gpr_seq_kernels(Xtrain, ytrain, k, sigma2, Xtest, ytest)
    cv_folds[fold,] <- t(gpr_seq_out)

  }
  cvm <- colMeans(cv_folds)
  return (list(cmv = cmv, k=k))
}
