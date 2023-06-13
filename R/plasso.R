#' Lasso and Post-Lasso
#' 
#' @description
#' \emph{plasso()} implicitly estimates a Lasso model using the \code{\link[glmnet]{glmnet}} package
#' and additionally estimates coefficient paths for a subsequent Post-Lasso model.
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y Vector of outcomes
#' @param w Vector of weights
#' @param ... Pass \code{\link[glmnet]{glmnet}} options
#' 
#' @import glmnet
#' @import Matrix
#' @importFrom stats coef predict var approx
#' @importFrom graphics axis matplot text
#' @importFrom methods as
#'
#' @return List including glmnet object and names of selected variables at cross-validated minima for Lasso and Post-Lasso.
#'
#' @export
#'
plasso = function(x,y,
                  w=NULL,
                  ...) {
  
  this.call = match.call()
  # Handle potentially provided sample weights, otherwise create weight vector of ones
  w = handle_weights(w,nrow(x))
  # Create variable names if not provided
  if ( is.null( colnames(x) ) ) colnames(x) = sprintf("var%s",seq(1:ncol(x)))
  
  # Lasso with full estimation sample
  lasso_full = glmnet(x,y,weights=as.vector(w),family="gaussian",...)
  coef_lasso_full = coef(lasso_full)                   # Save coefficients to extract later the ones at the CV minima
  lambda = lasso_full$lambda                           # Save grid to use the same in cross-validation
  x = add_intercept(x)
  
  l = length(lambda)
  coef_plasso_full = matrix(NA,
                            nrow=dim(coef_lasso_full)[1],ncol=dim(coef_lasso_full)[2],
                            dimnames=list(rownames(coef_lasso_full),colnames(coef_lasso_full)))
  
  for (i in 1:l) {
    
    nm_act = names(coef_lasso_full[,i])[which(coef_lasso_full[,i] != 0)]
    coef_plasso_full[,i] = fit_betas(x,y,w,nm_act,coef_lasso_full[,i])
    
  }
  coef_plasso_full = as(coef_plasso_full, "dgCMatrix")
  
  output = list("call"=this.call,
                "lasso_full"=lasso_full,
                "beta_plasso"=coef_plasso_full,
                "x"=x[,-1],"y"=y)
  
  class(output) = "plasso"
  return(output)
}


#' Plot coefficient paths
#' 
#' @description
#' Plot coefficient paths of (Post-) Lasso model.
#' 
#' @param x \code{\link[plasso]{plasso}} object
#' @param ... Pass generic \code{\link[base]{plot}} options
#' @param lasso If set as True, coefficient paths for Lasso instead of Post-Lasso is plotted. Default is False.
#' @param xvar What is on the X-axis
#' "norm" plots against the L1-norm of the coefficients,
#' "lambda" against the log-lambda sequence,
#' and "dev" against the percent deviance explained.
#' @param label If TRUE, label the curves with variable sequence numbers
#'
#' @export
#'
plot.plasso = function(x,..., lasso=FALSE,xvar=c("norm","lambda","dev"),label=FALSE) {
  
  
  nzC = utils::getFromNamespace("nonzeroCoef", "glmnet")
  
  if (!lasso) {
    
    lambda = x$lasso_full$lambda
    df = x$lasso_full$df
    dev = x$lasso_full$dev.ratio
    
    beta = x$beta_plasso[-1,]
    
    which = nzC(beta)
    nwhich = length(which)
    
    switch(nwhich+1,#we add one to make switch work
           "0"={
             warning("No plot produced since all coefficients zero")
             return()
           },
           "1"=warning("1 or less nonzero coefficients; glmnet plot is not meaningful")
    )
    beta = as.matrix(beta[which,,drop=FALSE])
    xvar = match.arg(xvar)
    switch(xvar,
           "norm"={
             index=apply(abs(beta),2,sum)
             iname="L1 Norm"
             approx.f=1
           },
           "lambda"={
             index=log(lambda)
             iname="Log Lambda"
             approx.f=0
           },
           "dev"= {
             index=dev
             iname="Fraction Deviance Explained"
             approx.f=1
           }
    )

    matplot(index,t(beta),lty=1,type="l",...)
    
    atdf = pretty(index)
    prettydf = approx(x=index,y=df,xout=atdf,rule=2,method="constant",f=approx.f)$y
    
    axis(3,at=atdf,labels=prettydf,tcl=NA)
    if(label){
      nnz = length(which)
      xpos = max(index)
      pos = 4
      if(xvar == "lambda"){
        xpos = min(index)
        pos = 2
      }
      xpos = rep(xpos,nnz)
      ypos = beta[,ncol(beta)]
      text(xpos,ypos,paste(which),cex=.5,pos=pos)
    }
    
  } else if (lasso) {
    
    plot(cv.plasso$lasso_full,xvar=xvar,label=label,...)
    
  }
  
}

