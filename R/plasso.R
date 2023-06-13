#' Lasso and Post-Lasso
#' 
#' @description
#' \emph{plasso()} implicitly estimates a Lasso model using the \code{\link{glmnet}} package
#' and additionally estimates coefficient paths for a subseqent Post-Lasso model.
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y Vector of outcomes
#' @param w Vector of weights
#' @param ... Pass \code{\link{glmnet}} options
#' 
#' @import glmnet
#' @importFrom stats coef predict var
#'
#' @return List including glmnet object and names of selected variables at cross-validated minima for Lasso and Post-Lasso.
#'
#' @export
#'
plasso = function(x,y,
                  w=NULL,
                  ...) {
  print(1)
}


#' Plot coefficient paths
#' 
#' @description
#' Plot coefficient paths of (Post-) Lasso model.
#' 
#' @param plasso \code{\link{plasso}} object
#' @param lasso If set as True, coefficient paths for Lasso instead of Post-Lasso is plotted. Default is False.
#'
#' @export
#'
plot.plasso = function(plasso,lasso=FALSE) {
  print(1)
  
}

