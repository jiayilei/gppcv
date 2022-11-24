source(utilities.R)
#' Title
#'
#' @param X training input data
#' @param y training response data
#' @param theta_init initlized parameter
#' @param iteration number of iteration for cross validation process, default to be 50
#' @param mode type of esimator used in the cv process, MAP is Maximum a Posteriori, MMSE is Minimum Mean Square Error
#'
#' @return best hyper parameter
#' @export
#'
#' @examples todo
pcv <- function(X, y, theta_init, iteration = 50, mode = c("MAP", "MMSE")){
  for (i in 1 : iteration){
    # split data set into half

    # obtain approximate mean function

    # find corresponding estimator theta based on the chosen mode

    # repeart procedure for the remaining half of the data

    # average thetas.

  }
  return (theta)

}
