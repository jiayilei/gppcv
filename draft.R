source("R/gp.R")
# test package
n <- 100
p <- 10
n_kernels <- 3

X <- matrix(rnorm(n * p, 0, 0.3), n, p)
y <- matrix(rnorm(n), n, 1)
num_folds <- 3
sigma2 <- 10


k <- c(se_kernel(3), se_kernel(3))
gpr_cv(X, y, k, sigma2, num_folds, NULL)


# pass in one kernel
k <- se_kernel(3)
gpr_cv(X, y, k, sigma2, num_folds, NULL)


#
library(standardGP)
library(matrixcalc)
# intent inputs
n <- 100 # number of inputs
p <- 10 # number of features
nt <- 50

set.seed(123)
X <- matrix(rnorm(n * p, 0, 0.3), n, p)
y <- matrix(rnorm(n), n, 1)
Xt <- matrix(rnorm(nt * p, 0, 0.3), nt, p)
yt <- matrix(rnorm(nt), nt, 1)

sigma2 <- 10
# first kernel, r is the Euclidean distance between two data points
k1 <- function(r) {
  return(exp(-r / 1.5))
}
# second kernel, r is the Euclidean distance between two data points
k2 <- function(r) {
  return(exp(-0.5 * (r / 1000)^2))
}
k <- c(k1, k2) # list of kernels

# function output

out <- gpr_seq_kernels(X, y, k, sigma2, Xt, yt)



pick_kernel(list(1), "se")
