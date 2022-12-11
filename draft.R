source("R/gp.R")
# test package
n = 100
p = 10
n_kernels = 3

X <- matrix(rnorm(n*p, 0, 0.3), n, p)
y = matrix(rnorm(n), n, 1)
num_folds=3
sigma2 = 10


k = c(se_kernel(3), se_kernel(3))
gpr_cv(X, y, k, sigma2, num_folds, NULL)


# pass in one kernel
k = se_kernel(3)
gpr_cv(X, y, k, sigma2, num_folds, NULL)

