#' Majority prediction for KTSP
#' 
#' Default prediction method for k-tsp.
#' @param x a matrix of predictions. Dimensions n_samples by k
#' @return vector of predictions
#' @family combinefun
#' @export 
majorityPredict <- function(x) {
    mean(x) < 0.5
}