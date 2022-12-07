
n = 10
p = 5
np = 5
nt = 3

X <- matrix(rnorm(n * p), n, p)
y <- matrix(rnorm(np), np)
Xt <- matrix(rnorm(nt * p), nt, p)

k = 1
sigma2 = 1
# check dimension
gpr(X, y, k, sigma2, Xt) # should stop


# test passing a function to a function
toy_func <- function(x,y){
  return (x + y)
}

use_func <- function(x, y, func){
  return (func(x,y))
}

x = 1
y = 2
use_func(x,y,toy_func)


# create a function that returns a function
A_func <- function(x){return (2*x)}

B_func <- function(x){return(x^2)}

method <- function(method){
  if (method == 'A'){
    return (function(x){
      return (2*x)
    })
  }
  if (method == 'B'){
    return (function(x){
      return (x^2)
    })
  }
}

method1 <-function(method){
  if (method == 'A'){
    return (A_func)
  }
  if (method == 'B'){
    return (B_func)
  }
}


test_fuc <- function(x, func){
  return (func(x))
}

K = method1('A')
x = 3
test_fuc(x, K)

# loop over a list of function
x = 3
fun_list <- c(A_func, B_func)
for (i in 1:length(fun_list)){
  temp = fun_list[[i]]
  print(temp(x))
}
