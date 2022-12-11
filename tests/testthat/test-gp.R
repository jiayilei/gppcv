library('tinytest')
# pre processing

# intent inputs
n = 100 # number of inputs
p = 10 # number of features
nt = 50 # number of test inputs
n_kernels = 3 # number of kernels

X <- matrix(rnorm(n*p, 0, 0.3), n, p)
Xt <- matrix(rnorm(nt*p, 0), nt, p)
y <- matrix(rnorm(n), n, 1)
yt <- matrix(rnorm(nt), nt, 1)

num_folds=3
fold_ids <- sample((1:n) %% num_folds + 1, n)
sigma2 = 10
k = c(se_kernel(3), se_kernel(3))


test_that("test covariance_mats", {
  # intent inputs
  n = 100 # number of inputs
  p = 10 # number of features
  nt = 50 # number of test inputs
  n_kernels = 3 # number of kernels

  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  Xt <- matrix(rnorm(nt*p, 0), nt, p)
  y <- matrix(rnorm(n), n, 1)
  yt <- matrix(rnorm(nt), nt, 1)

  num_folds=3
  fold_ids <- sample((1:n) %% num_folds + 1, n)
  sigma2 = 10
  k = c(se_kernel(3), se_kernel(3))

  kt <- k[[1]]
  co_out <- covariance_mats(X, Xt, kt)
  expect_equal(dim(co_out$K)[1], n)
  expect_equal(dim(co_out$K)[2], n)
  expect_equal(dim(co_out$Ks)[1], nt)
  expect_equal(dim(co_out$Ks)[2], p)
  expect_equal(co_out$K[1,2], se_kernel(2)(get_r(X[1,], X[2,])))
  expect_equal(co_out$K[3,2], se_kernel(2)(get_r(X[2,], X[3,])))
  expect_equal(co_out$ks[4,2], se_kernel(2)(get_r(Xt[4,], X[2,])))

  Xt_e <- matrix(rnorm(nt* 5, 0), nt, 5)
  expect_error(covariance_mats(X, Xt_e, k), "Dimensions of X and Xt do not match!")
})

test_that("test exp_kernel", {
  expect_equal(exp_kernel(2)(3), exp(-3/2))
  expect_error(exp_kernel(-1), "l needs to be positive")
  expect_error(exp_kernel(-1)(3), "l needs to be positive")
})

test_that("test get_r", {
  # intent inputs
  n = 100 # number of inputs
  p = 10 # number of features
  nt = 50 # number of test inputs
  n_kernels = 3 # number of kernels

  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  Xt <- matrix(rnorm(nt*p, 0), nt, p)
  y <- matrix(rnorm(n), n, 1)
  yt <- matrix(rnorm(nt), nt, 1)

  num_folds=3
  fold_ids <- sample((1:n) %% num_folds + 1, n)
  sigma2 = 10
  k = c(se_kernel(3), se_kernel(3))

  x1 = X[1,]
  x2 = X[2,]
  expect_equal(get_r(x1, x2), sqrt(sum((x1 - x2)^2)))
  x2 = c(3,4)
  expect_error(get_r(x1, x2), "Dimensions of x1 and x2 do not match!")
})


test_that("test gpr_cv", {
  # intent inputs
  n = 100 # number of inputs
  p = 10 # number of features
  nt = 50 # number of test inputs
  n_kernels = 3 # number of kernels

  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  Xt <- matrix(rnorm(nt*p, 0), nt, p)
  y <- matrix(rnorm(n), n, 1)
  yt <- matrix(rnorm(nt), nt, 1)

  num_folds=3
  fold_ids <- sample((1:n) %% num_folds + 1, n)
  sigma2 = 10
  k = c(se_kernel(3), se_kernel(3))

  # function output
  out <- gpr_cv(X, y, k, sigma2, num_folds, fold_ids)

  # manual calculation
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

  # check and compare
  expect_equal(out$cvm, cvm)

  ## test dimension check
  y_e <- matrix(rnorm(2*n), 2*n, 1)
  expect_error(gpr_cv(X, y_e, k, sigma2, num_folds, fold_ids), "Dimension of X and y does not match!")

})

test_that("test gpr_seq_kernels", {
  # intent inputs
  n = 100 # number of inputs
  p = 10 # number of features
  nt = 50 # number of test inputs
  n_kernels = 3 # number of kernels

  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  Xt <- matrix(rnorm(nt*p, 0), nt, p)
  y <- matrix(rnorm(n), n, 1)
  yt <- matrix(rnorm(nt), nt, 1)

  num_folds=3
  fold_ids <- sample((1:n) %% num_folds + 1, n)
  sigma2 = 10
  k = c(se_kernel(3), se_kernel(3))


  # inputs used to test dimension checks
  Xte <- matrix(rnorm( nt *5, 0, 0.3), nt, 5)
  ye <- matrix(rnorm(nt), nt, 1)
  yte <- matrix(rnorm(0.8*nt), 0.8*nt, 1)
  sigma2e <- -1

  # test dimension checks
  expect_error(gpr_seq_kernels(X, y, k, sigma2, Xte, yt), "Dimensions of X and Xt do not match!")
  expect_error(gpr_seq_kernels(X, ye, k, sigma2, Xt, yt), "Dimensions of X and y do not match!")
  expect_error(gpr_seq_kernels(X, y, k, sigma2, Xt, yte), "Dimensions of yt and Xt do not match!")
  expect_error(gpr_seq_kernels(X, y, k, sigma2e, Xt, yt), "Sigma2 needs to be positive!")

  # output from the function
  out <- gpr_seq_kernels(X, y, k, sigma2, Xt, yt)

  ##  manual calcuations
  # standardize inputs
  standardized_out <- standardize(X, y, Xt, yt)
  X = standardized_out$X
  y = standardized_out$y
  Xt = standardized_out$Xt
  yt = standardized_out$yt

  # standardized gaussian process regression for each k
  mse <- matrix(rep(0, length(k)))
  for (i in 1 : length(k)){
    if (length(k) == 1) {k_single = k} else {k_single = k[[i]]}
    cov_out <- covariance_mats(X, Xt, k_single)
    K = cov_out$K
    ks = cov_out$ks
    gpr_out <- gpr_standardized(X, y, k_single, sigma2, Xt, yt, K, ks)
    # evaluate mse for each kernels
    mse[i] = sum((gpr_out$fs - yt)^2)/nt
  }

  expect_equal(out$mse, mse)

})



test_that("test gpr_standardized", {
  # intent inputs
  n = 100 # number of inputs
  p = 10 # number of features
  nt = 50 # number of test inputs
  n_kernels = 3 # number of kernels

  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  Xt <- matrix(rnorm(nt*p, 0), nt, p)
  y <- matrix(rnorm(n), n, 1)
  yt <- matrix(rnorm(nt), nt, 1)

  num_folds=3
  fold_ids <- sample((1:n) %% num_folds + 1, n)
  sigma2 = 10
  k = se_kernel(3)

  cov_out <- covariance_mats(X, Xt, k)
  K = cov_out$K
  ks = cov_out$ks

  # inputs used to test dimension checks
  Ke <- matrix(rep(1, n*n),n, n)
  sigma2e <- 1e-10

  # test dimension checks
  ks = matrix(rep(1, nt*n),nt,n)
  expect_error(gpr_standardized(X, y, k, sigma2e, Xt, yt, Ke, ks), "Can not perform cholesky decomposition. Increase sigma2.")

  # output from the function
  out <- gpr_standardized(X, y, k, sigma2, Xt, yt, K, ks)

  ##  manual calcuations
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

  ## check
  expect_equal(out$fs, fs)
  expect_equal(out$Vfs, Vfs)
  expect_equal(out$logp, logp)

})



test_that("test matern_kernel", {
  expect_equal(matern_kernel(2,3)(4), 1 / gamma(3) / 2^(3-1) * (sqrt(2*3)/ 2 * 4)^ 3 * besselK(sqrt(2*3) * 4 / 2, nu = 3))
  expect_error(matern_kernel(-2,3)(4), "l needs to be positive")
  expect_error(matern_kernel(2,-3)(4), "v needs to be positive")
})
