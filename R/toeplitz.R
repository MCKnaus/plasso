#' Simulated 'Toeplitz' Data
#' 
#' @description
#' Some simulated data with a causal relationship between features X and target
#' y. The first 10 of the 25 features have a decreasing causal effect on y
#' while the remaining 15 do not have causal effect.
#' The variables in X follow a normal distribution with mean zero while
#' the covariance matrix follows a 'Toeplitz' matrix.
#'
#' @docType data
#' 
#' @usage data(toeplitz)
#'
#' @keywords datasets
#' 
#' @export
#' 
#' @examples
#' # load toeplitz data
#' data(toeplitz)
#' # extract target and features from data
#' y = as.matrix(toeplitz[,1])
#' X = toeplitz[,-1]
#' # fit cv.plasso to the data
#' \donttest{p.cv = plasso::cv.plasso(X,y)}
#'
"toeplitz"