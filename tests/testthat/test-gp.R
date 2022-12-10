library('tinytest')
test_that("test covariance_mats", {
  n = 300
  p = 10
  nt = 100
  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  Xt <- matrix(rnorm(nt*p, 0), nt, p)
  k = se_kernel(2)
  co_out <- covariance_mats(X, Xt, k)
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
  n = 300
  p = 10
  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  x1 = X[1,]
  x2 = X[2,]
  expect_equal(get_r(x1, x2), sqrt(sum((x1 - x2)^2)))
  x2 = c(3,4)
  expect_error(get_r(x1, x2), "Dimensions of x1 and x2 do not match!")
})


test_that("test gpr_cv", {
  n = 50 # number of inputs
  p = 10 # number of features
  n_kernels = 3 # number of kernels
  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  y <- matrix(rnorm(n), n, 1)
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
  # intended inputs
  n = 100
  p = 10
  nt = 50
  X <- matrix(rnorm(n*p, 0, 0.3), n, p)
  Xt <- matrix(rnorm(nt*p, 0), nt, p)
  y <- matrix(rnorm(n), n, 1)
  yt <- matrix(rnorm(nt), nt, 1)
  k = c(se_kernel(3), se_kernel(3))
  sigma2 <- 5

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

  # get K and Ks, covariance matrix

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

