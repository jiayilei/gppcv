
<!-- README.md is generated from README.Rmd. Please edit that file -->

# standardGP

<!-- badges: start -->

[![R-CMD-check](https://github.com/jiayilei/Project_Rpackage/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jiayilei/Project_Rpackage/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/jiayilei/Project_Rpackage/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jiayilei/Project_Rpackage?branch=main)
<!-- badges: end -->

## Description:

This package is based on predictions and log marginal likelihood for
Gaussian process regression (GP regression) algorithm proposed by C.
Rasmussen and C. William\[1\]. This package is intended to use the GP
regression algoithm to make prediciton based on gaussian-distribuated
mutilvariate data. User can pass in a sequence of different kernels and
compare the prediction accuracy across different kernels.

## Installation

You can install the development version of standardGP from
[GitHub](https://github.com/jiayilei/gppcv) with:

``` r
# install.packages("devtools")
devtools::install_github("jiayilei/gppcv")
```

Then, use

    library(standardGP)

## Examples

Below is a basic example which shows you how to perform cross validation
with a sequence of kernels by using `gpr_cv`:

``` r
library(standardGP)
library(matrixcalc)
# intent inputs
n <- 100 # number of inputs
p <- 10 # number of features
n_kernels <- 3 # number of kernels


set.seed(123)
X <- matrix(rnorm(n * p, 0, 0.3), n, p)
y <- matrix(rnorm(n), n, 1)

num_folds <- 3
fold_ids <- sample((1:n) %% num_folds + 1, n)
sigma2 <- 10
# first kernel, r is the Euclidena distance between two data points
k1 <- function(r) {
  return(exp(-r / 1.5))
}
# second kernel, r is the Euclidena distance between two data points
k2 <- function(r) {
  return(exp(-0.5 * (r / 1000)^2))
}
k <- c(k1, k2) # list of kernels

# function output
out <- gpr_cv(X, y, k, sigma2, num_folds, fold_ids)
out
#> $cvm
#> [1] 1.065387 1.073015
#> 
#> $k
#> $k[[1]]
#> function(r) {
#>   return(exp(-r / 1.5))
#> }
#> <bytecode: 0x0000029dbb7ecc78>
#> 
#> $k[[2]]
#> function(r) {
#>   return(exp(-0.5 * (r / 1000)^2))
#> }
#> <bytecode: 0x0000029dbbb10e98>
```

If you are interested in performing Gaussian process regression by using
multiple kernels, below is an example of how to use `gpr_seq_kernels` to
do the prediction. Note that kernels passed into this function need to
be using the Euclidean distance between two data points as parameters.

``` r
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
```

There are some predefined functions in the package that you can choose
from by using `pick_kernel`:

``` r
pick_kernel(list(1), "se")
#> function(r){
#>     return (exp(-0.5 * (r/l)^2))
#>   }
#> <bytecode: 0x0000029dbefcf788>
#> <environment: 0x0000029dbefda6b8>
pick_kernel(list(2, 3), "m")
#> function (r){
#>     left <-  1 / gamma(v) / 2^(v-1)
#>     mid <- (sqrt(2*v)/ l * r)^ v
#>     right <- besselK(sqrt(2*v) * r / l, nu = v)
#>     return (left * mid * right)
#>   }
#> <bytecode: 0x0000029dbf08a438>
#> <environment: 0x0000029dbf097f30>
```

## References

\[1\] Rasmussen, Carl Edward, and Christopher K. I. Williams. Gaussian
Processes for Machine Learning. MIT Press, 2005.
