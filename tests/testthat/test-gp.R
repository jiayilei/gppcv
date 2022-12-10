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


