#
#' Gaussian Process Regression with Standardized Inputs
#' @param X inputs
#' @param y targets
#' @param k covariance function
#' @param sigma2 noise level
#' @param Xt test inputs
#' @param yt test targets
#' @param K covariance matrix between training data
#' @param ks covariance matrix between testing data and training data
#'
#' @return
#'        fs mean function
#'        Vfs variance
#'        logp log marginal likelihood
#' @export
#'
#' @examples
#' n = 10
#' p = 3
#' nt = 4
#' X <- matrix(rnorm(n*p, 0, 0.3), n, p)
#' y = matrix(rnorm(n), n, 1)
#' Xt <- matrix(rnorm(nt*p, 0, 0.3), nt, p)
#' yt = matrix(rnorm(nt), nt, 1)
#' num_folds=3
#' sigma2 = 10
#'
#'
#' k <- function(r){
#'   return (exp(-0.5 * (r/4)^2))
#' }
#'
#' # covariance matrix: covariance evaluated at all pairs of training point
#'  K = matrix(rep(0,n*n), n,n)
#'  # covariance vector: covariance between test point and the n training points
#'  ks = matrix(rep(0,nt * n), nt, n)

#'  for (i in (1:n)){
#'    for (j in (1:n)){
#'      r = get_r(X[i,], X[j,])
#'      K[i,j] = k(r)
#'    }
#'    for (m in (1: nt)){
#'      r = get_r(X[i,], Xt[m,])
#'      ks[m,i] = k(r)
#'    }
#'  }
#'  tK = t(K)
#'  diag(tK) <- 0
#'  K = tK + K
#'
#'
#' gpr_standardized(X, y, k, sigma2, Xt, yt, K, ks)
#'
gpr_standardized <- function(X, y, k, sigma2, Xt, yt, K, ks){
  # compatibility check
  n = dim(X)[1]
  p = dim(X)[2]
  nt = dim(Xt)[1]

  # cholesky factorization
  # positive definit check
  if (!(is.positive.definite(K + sigma2 * diag(n), tol=1e-10))){
    stop("Can not perform cholesky decomposition. Increase sigma2.")
  }
  L = chol(K + sigma2 * diag(n))

  # find predictive mean
  alpha = solve(t(L)) %*% (solve(L) %*% y)
  fs = ks %*% alpha

  # find predictive variance
  v = solve(L) %*% t(ks)
  r = matrix(rep(0, nt), nt)
  for (i in 1:nt){
    r[i] = get_r(Xt[i,], Xt[i,])
  }
  Vfs = k(r) - colSums(v^2)

  # calculate log marginal likelihood
  logp = -1/2 * crossprod(y, alpha) - sum(log(diag(L))) - n/2 * log(2*pi)

  return (list(fs = fs, Vfs = Vfs, logp = logp))
}


#' Standardize Input
#' Standardize X, y, Xt, and yt based on mean and variance of X and y
#' @param X original training inputs
#' @param y original training targets
#' @param Xt original testing inputs
#' @param yt original testing targets
#'
#' @return
#'        X scaled training inputs
#'        y scaled training targets
#'        Xt scaled testing inputs
#'        yt scaled testing targets
#' @export
#'
#' @examples
#' n = 10
#' p = 3
#' nt = 4
#' X <- matrix(rnorm(n*p, 0, 0.3), n, p)
#' y = matrix(rnorm(n), n, 1)
#' Xt <- matrix(rnorm(nt*p, 0, 0.3), nt, p)
#' yt = matrix(rnorm(nt), nt, 1)
#' standardize(X, y, Xt, yt)
#'
#'
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


#
#' Returns Euclidean distance between two input data
#'
#' @param x1 first data
#' @param x2 second data
#'
#' @return r euclidean distance between first data and the second data
#' @export
#'
#' @examples
#' x1 = matrix(rnorm(9))
#' x2 = matrix(rnorm(9))
#' get_r(x1, x2)
get_r <- function(x1, x2){
  if (length(x1) != length(x2)){
    stop("Dimensions of x1 and x2 do not match!")
  }
  return (sqrt(sum((x1 - x2)^2)))
}


#
#' Square Exponential Kernel
#'
#' @param l length scale parameter of squared exponential kernel
#' @param r euclidean distance between two data points
#'
#' @return covariance function with fixed hyperperameter l
#' @export
#'
#' @examples
#' l = 3
#' se_kernel(3)
se_kernel <- function(l, r = 1){
  if (l <=0){stop("l needs to be positive")}
  fun <- function(r){
    return (exp(-0.5 * (r/l)^2))
  }
  return (fun)
}


#' Matern kernel
#'
#' @param l length-scale
#' @param v smoothness
#' @param r euclidean distance between two data points
#'
#' @return covariance function with fixed parameter l and v
#' @export
#'
#' @examples
#' matern_kernel(3, 2.5)
#' matern_kernel(1, 1.5)
matern_kernel <- function(l, v, r = 1){
  # input check
  if (l <=0){stop("l needs to be positive")}
  if (v <=0){stop("v needs to be positive")}
  fun <- function (r){
    left <-  1 / gamma(v) / 2^(v-1)
    mid <- (sqrt(2*v)/ l * r)^ v
    right <- besselK(sqrt(2*v) * r / l, nu = v)
    return (left * mid * right)
  }
  return (fun)
}


#' Exponential kernel
#'
#' @param l length scale
#' @param r euclidean distance between two data points
#'
#' @return covariance function with fixed parameter l
#' @export
#'
#' @examples
#' exp_kernel(3)
#' exp_kernel(3.5)
exp_kernel <- function(l, r=1){
  # input check
  if (l <=0){stop("l needs to be positive")}
  fun <- function(r){
    return (exp(-r/l))
  }
  return (fun)
}


#
#' pick the kernel
#' return kernel functions with hyperparemeters
#' @param method a string indicates what function class to pick
#' @param para_list a list of parameter that will be pass on to the kernel function
#'
#' @return return one of the kernel functions with input parameter list
#' @export
#'
#' @examples
#' pick_kernel(list(1,2), 'm')
#' pick_kernel(list(3), 'se')
pick_kernel <- function (para_list, method = c('se', 'm', 'exp')){
  method = match.arg(method)
  if (method == 'se'){
    if (length(para_list) != 1){
      stop("squared exponential kernel needs 1 parameter")
    }
    return (do.call(se_kernel, para_list))
  } else if (method == 'm'){
    if (length(para_list) != 2){
      stop("matern kernel needs 2 parameters")
      }
    return (do.call(matern_kernel, para_list))
  } else if (method == 'exp'){
    if (length(para_list) != 1){
      stop("exponential kernel needs 1 parameter")
    }
    return (do.call(exp_kernel, para_list))
  }

}

#
#' Covariance Matrix
#' gets covariance_matrix of K and Ks
#'
#' @param X training inputs
#' @param Xt testing inputs
#' @param k covariance function
#'
#' @return
#'     K covariance matrix of training input
#'     Ks covariance matrix of trainign data and test data
#' @export
#'
#' @examples
#' k <- function(r){
#'   return (exp(-0.5 * (r/4)^2))
#' }
#' n = 10
#' p = 3
#' nt = 4
#' X <- matrix(rnorm(n*p, 0, 0.3), n, p)
#' Xt <- matrix(rnorm(nt*p, 0, 0.3), nt, p)
#'
covariance_mats <- function(X, Xt, k){

  n = dim(X)[1]
  p = dim(X)[2]
  nt = dim(Xt)[1]

  if (p != dim(Xt)[2]){
    stop("Dimensions of X and Xt do not match!")
  }

  # covariance matrix: covariance evaluated at all pairs of training point
  K = matrix(rep(0,n*n), n,n)
  # covariance vector: covariance between test point and the n training points and variance of test data points
  ks = matrix(rep(0, nt * n), nt, n)

  for (i in (1:n)){
    # covariance between training data
    for (j in (i:n)){
      r = get_r(X[i,], X[j,])
      K[i,j] = k(r)
    }
    # covariance between test data and training data
    for (m in (1: nt)){
      r = get_r(X[i,], Xt[m,])
      ks[m,i] = k(r)
    }
  }

  tK = t(K)
  diag(tK) <- 0
  K = tK + K

  return (list(K = K, ks = ks))
}


# sequence of k's to fit the model and compare different kernels performance
#' Gaussian Process Regression for a sequence of kernels
#'
#' @param X training inputs
#' @param y training targets
#' @param k list of covariance functions
#' @param sigma2 noise level
#' @param Xt testing inputs
#' @param yt testing targets
#'
#' @return
#'        k: list of covariance functions, same as input
#'        mse: mean square errors for each covariance function
#'        predictions: [nt x num_kernels] matrix containing predictions of each test inputs under each kernels
#' @export
#'
#' @examples
#' n = 10
#' p = 3
#' nt = 4
#' X <- matrix(rnorm(n*p, 0, 0.3), n, p)
#' y = matrix(rnorm(n), n, 1)
#' Xt <- matrix(rnorm(nt*p, 0, 0.3), nt, p)
#' yt = matrix(rnorm(nt), nt, 1)
#' num_folds=3
#' sigma2 = 10
#'
#'
#' k1 <- function(r){
#'   return (exp(-0.5 * (r/4)^2))
#' }
#'
#' k2 <- function(r){
#'   return (exp(-r/3.5))
#' }
#' k = c(k1, k2)
#' gpr_seq_kernels(X, y, k, sigma2, Xt, yt)
#'
gpr_seq_kernels <-function(X, y, k, sigma2, Xt, yt){
  # compatibility check
  n = dim(X)[1]
  p = dim(X)[2]
  nt = dim(Xt)[1]

  if (length(y) != n){
    stop("Dimensions of X and y do not match!")
  }
  if (dim(Xt)[2] != p){
    stop("Dimensions of X and Xt do not match!")
  }
  if (length(yt) != nt){
    stop("Dimensions of yt and Xt do not match!")
  }
  if (sigma2 < 0){
    stop("Sigma2 needs to be positive!")
  }

  # standardize inputs
  standardized_out <- standardize(X, y, Xt, yt)
  X = standardized_out$X
  y = standardized_out$y
  Xt = standardized_out$Xt
  yt = standardized_out$yt

  # standardized gaussian process regression for each k
  mse <- matrix(rep(0, length(k)))
  predictions <- matrix(rep(0, nt*length(k)), nt, length(k))

  for (i in 1 : length(k)){
    if (length(k) == 1) {k_single = k} else {k_single = k[[i]]}
    cov_out <- covariance_mats(X, Xt, k_single)
    K = cov_out$K
    ks = cov_out$ks
    gpr_out <- gpr_standardized(X, y, k_single, sigma2, Xt, yt, K, ks)
    # evaluate mse for each kernels
    mse[i] = sum((gpr_out$fs - yt)^2)/nt
    predictions[,i] = gpr_out$fs
  }
  return (list(k=k, mse = mse, predictions = predictions))

}

#
#' Gaussian Process Regression Cross Validation
#'
#' @param X training inputs
#' @param y training targets
#' @param k list of kernels
#' @param sigma2 noise level
#' @param num_folds number of folds
#' @param fold_ids list of ids to indicate folds
#'
#' @return
#'        cvm average mean square errors for each kernels
#'        k original kernels
#' @export
#'
#' @examples
#' n = 100
#' p = 3
#' nt = 4
#' X <- matrix(rnorm(n*p, 0, 0.3), n, p)
#' y = matrix(rnorm(n), n, 1)
#' sigma2 = 10
#' k1 <- function(r){
#'   return (exp(-0.5 * (r/4)^2))
#' }
#'
#' k2 <- function(r){
#'   return (exp(-r/3.5))
#' }
#' k = c(k1, k2)
#' num_folds = 5
#' gpr_cv(X, y, k, sigma2, num_folds)
gpr_cv <- function(X, y, k, sigma2, num_folds, fold_ids = NULL){
  n = dim(X)[1]
  if (n != dim(matrix(y))[1]){
    stop("Dimension of X and y does not match!")
  }
  # get a list of fold_ids if it is not given
  if (is.null(fold_ids) && !(is.null(num_folds))) {
    fold_ids <- sample((1:n) %% num_folds + 1, n)
  }else {
    num_folds = max(fold_ids)
  }

  nkernels <- length(k)
  cv_folds <- matrix(rep(0, num_folds * nkernels), num_folds, nkernels)
  for (fold in 1:num_folds) {
    # split data into training set and testing set
    Xtrain <- X[fold_ids != fold, ]
    ytrain <- y[fold_ids != fold]

    Xtest <- X[fold_ids == fold, ]
    ytest <- y[fold_ids == fold]

    # train models with sequence of kernels
    gpr_seq_out <- gpr_seq_kernels(Xtrain, ytrain, k, sigma2, Xtest, ytest)
    # browser()
    cv_folds[fold,] <- t(gpr_seq_out$mse)

  }
  # take avg. mse of each kernels
  cvm <- colMeans(cv_folds)
  return (list(cvm = cvm, k=k))
}
