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


#' Print (Post-) Lasso model
#' 
#' @description
#' Printing main insights from (Post-) Lasso model.
#'
#' @param x \code{\link[plasso]{plasso}} object
#' @param ... Pass generic \code{\link[base]{print}} options
#' @param digits Integer, used for number formatting
#'
#' @return Prints glmnet like output.
#' 
#' @method print plasso
#'
#' @export
#'
print.plasso = function(x,...,digits=max(3, getOption("digits")-3)) {
  
  plasso = x
  cat("\nCall: ", deparse(plasso$call), "\n\n")
  out = data.frame(Df=plasso$lasso_full$df,
                   `%Dev`=round(plasso$lasso_full$dev.ratio*100, 2),
                   Lambda=signif(plasso$lasso_full$lambda, digits),
                   check.names=FALSE,row.names=seq(along=plasso$lasso_full$df))
  class(out) = c("anova",class(out))
  print(out)
}


#' Summary of (Post-) Lasso model
#' 
#' @description
#' Summary of (Post-) Lasso model.
#'
#' @param object \code{\link[plasso]{plasso}} object
#' @param ... Pass generic \code{\link[base]{summary}} summary options
#'
#' @return 'summaryDefault' object
#' 
#' @method summary plasso
#'
#' @export
#'
summary.plasso = function(object,...) {
  
  return(summary.default(object, ...))
  
}


#' Predict for (Post-) Lasso models
#' 
#' @description
#' Prediction for (Post-) Lasso models.
#'
#' @param object Fitted \code{\link[plasso]{plasso}} model object
#' @param ... Pass generic \code{\link[stats]{predict}} options
#' @param xnew Matrix of new values for x at which predictions are to be made. If no value is supplied, x from fitting procedure is used. This argument is not used for type="coefficients".
#' @param type Type of prediction required. "response" returns fitted values, "coefficients" returns beta estimates.
#' @param s If Null, prediction is done for all lambda values. If a value is provided, the closest lambda value of the plasso object is used.
#' @param weights Vector of weights. This argument is not used for default type="response".
#' 
#' @return Returns predictions of coefficients or fitted values for all lambda values or a particular one.
#'
#' @method predict plasso
#'
#' @export
#'
predict.plasso = function(object,
                          ...,
                          xnew=NULL,
                          type=c("response","coefficients"),
                          s=NULL,
                          weights=NULL
                          ) {
  
  plasso = object
  type = match.arg(type)
  
  if (is.null(xnew)) x = plasso$x else x = xnew
  if ( is.null( colnames(x) ) ) colnames(x) = sprintf("var%s",seq(1:ncol(x)))
  x = add_intercept(x)
  y = plasso$y
  w = handle_weights(w,nrow(x))
  
  if (is.null(s)) {
    
    l = length(plasso$lasso_full$lambda)
    
    if (type == "coefficients") {
      
      coef_lasso = matrix(NA,nrow=l,ncol=ncol(x),dimnames=list(1:l,colnames(x)))
      coef_plasso = matrix(NA,nrow=l,ncol=ncol(x),dimnames=list(1:l,colnames(x)))
      
      for (i in 1:l) {
        coef_lasso[i,] = coef(plasso$lasso_full)[,i]
        
        nm_act = names(coef(plasso$lasso_full)[,i])[which(coef(plasso$lasso_full)[,i] != 0)]
        coef_plasso[i,] = fit_betas(x,y,w,nm_act,coef(plasso$lasso_full)[,i])
      }
      colnames(coef_lasso) = colnames(x)
      colnames(coef_plasso) = colnames(x)
      
      return(list("lasso"=coef_lasso,"plasso"=coef_plasso))
      
    } else if (type == "response"){
      
      fit_lasso = matrix(NA,nrow=nrow(x),ncol=l,dimnames=list(1:nrow(x),1:l))
      fit_plasso = matrix(NA,nrow=nrow(x),ncol=l,dimnames=list(1:nrow(x),1:l))
      
      for (i in 1:l) {
        
        fit_lasso[,i] = x %*% coef(plasso$lasso_full)[,i]
        
        nm_act = names(coef(plasso$lasso_full)[,i])[which(coef(plasso$lasso_full)[,i] != 0)]
        coef_plasso = fit_betas(x,y,w,nm_act,coef(plasso$lasso_full)[,i])
        fit_plasso[,i] = x %*% coef_plasso 
      }
      
      return(list("lasso"=fit_lasso,"plasso"=fit_plasso))
    }
      
  } else if (is.numeric(s) && length(s) == 1){
    
    abs_diff <- abs(plasso$lasso_full$lambda - s)
    closest_index <- which.min(abs_diff)
    closest_lambda <- plasso$lasso_full$lambda[closest_index]
    
    coef_lasso = coef(plasso$lasso_full)[,closest_index]
    
    nm_act = names(coef_lasso)[which(coef_lasso != 0)]
    coef_plasso = fit_betas(x,y,w,nm_act,coef_lasso)
    
    if (type == "coefficients") {
      
      return(list("lasso"=coef_lasso,"plasso"=coef_plasso))
      
    } else if (type == "response"){
      
      fit_lasso = x %*% coef_lasso
      fit_plasso = x %*% coef_plasso 
      
      return(list("lasso"=fit_lasso,"plasso"=fit_plasso))
      
    }
  }
}


#' Extract coefficients from a plasso object
#' 
#' @description
#' Extract coefficients for both Lasso and Post-Lasso from a plasso object.
#' 
#' @param object \code{\link[plasso]{plasso}} object
#' @param ... Pass generic \code{\link[stats]{coef}} options
#' @param s If Null, coefficients are returned for all lambda values. If a value is provided, the closest lambda value of the plasso object is used.
#' @param weights Vector of weights for fitting coefficients.
#' 
#' @return List containing matrices or vector of coefficients for both Lasso and Post-Lasso.
#'
#' @method coef plasso
#'
#' @export 
#' 
coef.plasso = function(object,...,s=NULL,weights=NULL){
  return(predict(object,...,s=s,type="coefficients",weights=weights))
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
#' @method plot plasso
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

