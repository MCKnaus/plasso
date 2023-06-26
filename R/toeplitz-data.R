#' Simulated data based on toeplitz matrix
#'
#' Simulated data relying on a toeplitz like variance-covariance matrix of
#' features X
#'
#' @docType data
#'
#' @usage data(toeplitz)
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(toeplitz)
#' y = as.matrix(toeplitz[,1])
#' X = toeplitz[,-1]
#' \donttest{p = plasso::cv.plasso(X,y)}