
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

This is a basic example which shows you how to perform cross valdiation
with a sequence of kernels by using `gpr_cv`:

``` r
library(standardGP)
library(matrixcalc)
# intent inputs
n = 100 # number of inputs
p = 10 # number of features
n_kernels = 3 # number of kernels

set.seed(123)
X <- matrix(rnorm(n*p, 0, 0.3), n, p)
y <- matrix(rnorm(n), n, 1)

num_folds=3
fold_ids <- sample((1:n) %% num_folds + 1, n)
sigma2 = 10
# first kernel
k1  <- function(r){
    return (exp(-r/1.5))
}
# second kernel
k2 <- function(r){
    return (exp(-0.5 * (r/1000)^2))
}
k = c(k1, k2) # list of kernels

# function output
out <- gpr_cv(X, y, k, sigma2, num_folds, fold_ids)
out
#> $cvm
#> [1] 1.065387 1.073015
#> 
#> $k
#> $k[[1]]
#> function(r){
#>     return (exp(-r/1.5))
#> }
#> <bytecode: 0x000001f843839d60>
#> 
#> $k[[2]]
#> function(r){
#>     return (exp(-0.5 * (r/1000)^2))
#> }
#> <bytecode: 0x000001f843b67470>
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
