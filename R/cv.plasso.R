#' Cross-validated Lasso and Post-Lasso
#' 
#' @description
#' \emph{cv.plasso()} uses the \code{\link[glmnet]{glmnet}} package to estimate the coefficient paths and cross-validates least squares Lasso AND Post-Lasso.
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y Vector of outcomes
#' @param w Vector of weights
#' @param kf Number of folds in k-fold CV
#' @param parallel Set as TRUE for parallelized cross-validation. Default is FALSE.
#' @param ... Pass \code{\link[glmnet]{glmnet}} options
#' 
#' @import glmnet
#' @import Matrix
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import iterators
#' @importFrom stats coef predict var
#'
#' @return List object including base \code{\link[glmnet]{glmnet}} object and cross-validation results (incl. optimal Lambda values) for both Lasso and Post-Lasso model.
#'
#' @export
#'
cv.plasso = function(x,y,
                  w=NULL,
                  kf=10,
                  parallel=FALSE,
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
  
  # Get indicator for CV samples
  split = stats::runif(nrow(x))
  cvgroup = as.numeric(cut(split,stats::quantile(split,probs=seq(0,1,1/kf)),include.lowest=TRUE))  # Groups for k-fold CV
  list = 1:kf                                         # Needed later in the loop to get the appropriate samples
  
  if (!parallel) {
    
    # Initialize matrices for MSE of Lasso and post-lasso for each grid point and CV fold
    cv_MSE_lasso = matrix(nrow = kf,ncol = length(lambda))
    cv_MSE_plasso = matrix(nrow = kf,ncol = length(lambda))
    
    # Start loop for cross-validation of Lasso and Post-Lasso
    for (i in 1:kf) {
      CV = CV_core(x,y,w,cvgroup,list,i,lambda,...)
      
      # Extract MSE of Lasso and Post-Lasso
      cv_MSE_lasso[i,] = CV$MSE_lasso
      cv_MSE_plasso[i,] = CV$MSE_plasso
    }
    
  } else if (parallel)  {
    
    n_cores <- min(parallel::detectCores()-1,kf)
    
    cl = parallel::makeCluster(n_cores, methods=FALSE)
    doParallel::registerDoParallel(cl)
    
    para_res <- foreach::foreach(i = 1:kf, .combine="rbind", .inorder=FALSE) %dopar% {
                           
                           ## core CV program functionsxxx no. 9
                           CV = CV_core(x,y,w,cvgroup,list,i,lambda,...)
                           
                           ## Extract MSE of Lasso and Post-Lasso
                           cv_MSE_lasso = CV$MSE_lasso
                           cv_MSE_post_lasso = CV$MSE_plasso
                           
                           ## Matrices to be returned from cores
                           return(list(as.matrix(cv_MSE_lasso),as.matrix(cv_MSE_post_lasso)))
    }
    parallel::stopCluster(cl)
    
    para_res <- as.matrix(do.call(cbind,para_res))
    cv_MSE_lasso <- t(para_res[,1:kf])
    cv_MSE_plasso <- t(para_res[,(kf+1):(2*kf)])
    
  }
  
  # Calculate mean MSEs over all folds
  cv_MSE_lasso[is.na(cv_MSE_lasso)] = max(cv_MSE_lasso) # can happen if glmnet does not go over the full grid
  cv_MSE_plasso[is.na(cv_MSE_plasso)] = max(cv_MSE_plasso) # and/or Post-Lasso has not full rank
  mean_MSE_lasso = colMeans(cv_MSE_lasso)
  mean_MSE_plasso = colMeans(cv_MSE_plasso)
  
  # Get grid position of minimum MSE
  ind_min_l = which.min(mean_MSE_lasso)
  ind_min_pl = which.min(mean_MSE_plasso)
  # Get variable names at minima
  names_l = names(coef_lasso_full[,ind_min_l])[which(coef_lasso_full[,ind_min_l] != 0)]
  names_pl = names(coef_lasso_full[,ind_min_pl])[which(coef_lasso_full[,ind_min_pl] != 0)]
  
  # Get Lasso coefficients at minimum of Lasso
  coef_min_l = coef_lasso_full[,ind_min_l][which(coef_lasso_full[,ind_min_l] != 0)]
  # Get Lasso coefficients at minimum of Post-Lasso
  coef_min_pl = coef_lasso_full[,ind_min_pl][which(coef_lasso_full[,ind_min_pl] != 0)]
  
  # Return names and coefficients
  output = list("call"=this.call,
                "lasso_full"=lasso_full,"kf"=kf,
                "cv_MSE_lasso"=cv_MSE_lasso,"cv_MSE_plasso"=cv_MSE_plasso,
                "mean_MSE_lasso"=mean_MSE_lasso, "mean_MSE_plasso"=mean_MSE_plasso,
                "ind_min_l"=ind_min_l,"ind_min_pl"=ind_min_pl,
                "lambda_min_l"=lambda[ind_min_l],"lambda_min_pl"=lambda[ind_min_pl],
                "names_l"=names_l,"names_pl"=names_pl,
                "coef_min_l"=coef_min_l,"coef_min_pl"=coef_min_pl,
                "x"=x,"y"=y)
  
  class(output) = "cv.plasso"
  return(output)
}


#' Print cross-validated (Post-) Lasso model
#' 
#' @description
#' Printing main insights from cross-validated (Post-) Lasso model.
#'
#' @param x \code{\link[plasso]{cv.plasso}} object
#' @param ... Pass generic \code{\link[base]{print}} options
#' @param digits Integer, used for number formatting
#'
#' @return Prints cross-validated MSE for both Lasso and Post-Lasso model as well as the number of active variables.
#' 
#' @method print cv.plasso
#'
#' @export
#'
print.cv.plasso = function(x,...,digits=max(3, getOption("digits")-3)) {
  
  plasso = x
  cat("\nCall: ", deparse(plasso$call), "\n\n")
  out = data.frame(Df=plasso$lasso_full$df,
                   `%Dev`=round(plasso$lasso_full$dev.ratio*100, 2),
                   Lambda=signif(plasso$lasso_full$lambda, digits),
                   MSE_lasso=plasso$mean_MSE_lasso,
                   MSE_plasso=plasso$mean_MSE_plasso,
                   check.names=FALSE,row.names=seq(along=plasso$lasso_full$df))
  class(out) = c("anova",class(out))
  print(out)
}


#' Summary of cross-validated (Post-) Lasso model
#' 
#' @description
#' Summary of cross-validated (Post-) Lasso model.
#'
#' @param object \code{\link[plasso]{cv.plasso}} object
#' @param ... Pass generic \code{\link[base]{summary}} summary options
#' @param default TRUE for glmnet-like summary output, FALSE for more specific summary information
#'
#' @return For specific summary information: List object containing optimal lambda value and associated MSE for both cross-validated Lasso and Post-Lasso model.
#' 
#' @method summary cv.plasso
#'
#' @export
#'
summary.cv.plasso = function(object,...,default=TRUE) {
  
  plasso = object
  
  if (default) {
    
    return(summary.default(plasso, ...))
    
  } else {
  
    value = list(
      call = plasso$call,
      mean_MSE_lasso = plasso$mean_MSE_lasso,
      mean_MSE_plasso = plasso$mean_MSE_plasso,
      lasso_full = plasso$lasso_full,
      ind_min_lasso = plasso$ind_min_l,
      ind_min_plasso = plasso$ind_min_pl
    )
    class(value) = "summary.cv.plasso"
    return(value)

  }
}


#' Print summary of (Post-) Lasso model
#'
#' @description
#' Prints summary information of cv.plasso object
#'
#' @param x Summary of plasso object (either of class "summary.plasso' or "summaryDefault")
#' @param ... Pass generic R \code{\link[base]{print}} options
#' @param digits Integer, used for number formatting
#' 
#' @method print summary.cv.plasso
#'
#' @export
#' 
print.summary.cv.plasso = function(x,...,digits=max(3L, getOption("digits") - 3L)) {
  
  if (inherits(x,'summaryDefault')) {
    
    print.summaryDefault(x, digits=digits, ...)
    
  } else if (inherits(x,'summary.cv.plasso')){
    
    cat("\nCall:\n ", paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep = "")
    
    cat("Lasso:\n Minimum CV MSE Lasso: ",toString(signif(min(x$mean_MSE_lasso,na.rm=TRUE),digits)))
    cat("\n Lambda at minimum: ",toString(signif(x$lasso_full$lambda[x$ind_min_l],digits)))
    cat("\n Active variables at minimum: ",names(coef(x$lasso_full)[,x$ind_min_l])[which(coef(x$lasso_full)[,x$ind_min_l] != 0)])
    
    cat("\nPost-Lasso:\n Minimum CV MSE Post-Lasso: ",toString(signif(min(x$mean_MSE_plasso,na.rm=TRUE),digits)))
    cat("\n Lambda at minimum: ",toString(signif(x$lasso_full$lambda[x$ind_min_pl],digits)))
    cat("\n Active variables at minimum: ",names(coef(x$lasso_full)[,x$ind_min_pl])[which(coef(x$lasso_full)[,x$ind_min_pl] != 0)])
    
  }
  
}


#' Predict after cross-validated (Post-) Lasso
#' 
#' @description
#' Prediction for cross-validated (Post-) Lasso.
#'
#' @param object Fitted \code{\link[plasso]{cv.plasso}} model object
#' @param ... Pass generic \code{\link[stats]{predict}} options
#' @param xnew Matrix of new values for x at which predictions are to be made. If no value is supplied, x from fitting procedure is used. This argument is not used for type="coefficients".
#' @param type Type of prediction required. "response" returns fitted values, "coefficients" returns beta estimates.
#' @param s Determines whether prediction is done for all values of lambda ("all") or only for the optimal lambda ("optimal") according to the standard error-rule.
#' @param weights Vector of weights. This argument is not used for default type="response".
#' @param se_rule Only If equal to zero predictions from CV minimum (default). Negative values go in the direction of smaller
#' models (e.g. se_rule=-1 creates the standard 1SE rule), positive values go in the direction of larger models
#' (e.g. se_rule=1 creates the standard 1SE+ rule). This argument is not used for s="all".
#' 
#' @return Returns predictions of coefficients or fitted values for all lambda values or the optimal one.
#'
#' @method predict cv.plasso
#'
#' @export
#'
predict.cv.plasso = function(object,
                             ...,
                             xnew=NULL,
                             type=c("response","coefficients"),
                             s=c("optimal","all"),
                             weights=NULL,
                             se_rule=0) {
  
  plasso = object
  # Check if type and s is valid
  type = match.arg(type)
  s = match.arg(s)
  
  if (is.null(xnew)) x = plasso$x else x = xnew
  if ( is.null( colnames(x) ) ) colnames(x) = sprintf("var%s",seq(1:ncol(x)))
  x = add_intercept(x)
  y = plasso$y
  w = handle_weights(w,nrow(x))
  
  if (s == "optimal") {
    
    oneSE_lasso = sqrt(apply(plasso$cv_MSE_lasso,2,var)/plasso$kf)
    oneSE_plasso = sqrt(apply(plasso$cv_MSE_plasso,2,var)/plasso$kf)
    oneSE_lasso[is.na(oneSE_lasso)] = 0
    oneSE_plasso[is.na(oneSE_plasso)] = 0
    ind_Xse_l = find_Xse_ind(plasso$mean_MSE_lasso,plasso$ind_min_l,oneSE_lasso,se_rule)
    ind_Xse_pl = find_Xse_ind(plasso$mean_MSE_plasso,plasso$ind_min_pl,oneSE_plasso,se_rule)
    
    coef_lasso = coef(plasso$lasso_full)[,ind_Xse_l]
    
    nm_act = names(coef(plasso$lasso_full)[,ind_Xse_pl])[which(coef(plasso$lasso_full)[,ind_Xse_pl] != 0)]
    coef_plasso = fit_betas(x,y,w,nm_act,coef(plasso$lasso_full)[,ind_Xse_pl])
    
    if (type == "coefficients") {
      
      return(list("lasso"=coef_lasso,"plasso"=coef_plasso))
      
    } else if (type == "response"){
      
      # Fitted values for lasso
      fit_lasso = x %*% coef_lasso
      
      # Fitted values for post lasso
      fit_plasso = x %*% coef_plasso 
      
      return(list("lasso"=fit_lasso,"plasso"=fit_plasso))
      
    }
    
  } else if (s == "all"){
    
    l = length(plasso$lasso_full$lambda)
    
    coef_lasso = as.matrix(coef(plasso$lasso_full))
    coef_plasso = matrix(NA,nrow=ncol(x),ncol=l,dimnames=list(colnames(x),colnames(coef_lasso)))
    
    for (i in 1:l) {
  
      nm_act = names(coef_lasso[,i])[which(coef_lasso[,i] != 0)]
      coef_plasso[,i] = fit_betas(x,y,w,nm_act,coef_lasso[,i])
    }
    
    if (type == "coefficients") {
      
      return(list("lasso"=coef_lasso,"plasso"=coef_plasso))
      
    } else if (type == "response"){
      
      fit_lasso = matrix(NA,nrow=nrow(x),ncol=l,dimnames=list(1:nrow(x),colnames(coef_lasso)))
      fit_plasso = matrix(NA,nrow=nrow(x),ncol=l,dimnames=list(1:nrow(x),colnames(coef_lasso)))
      
      for (i in 1:l) {
        
        fit_lasso[,i] = x %*% coef_lasso[,i]
        fit_plasso[,i] = x %*% coef_plasso[,i]
      }
      
      return(list("lasso"=fit_lasso,"plasso"=fit_plasso))
      
    }
  }
}


#' Extract coefficients from a cv.plasso object
#' 
#' @description
#' Extract coefficients for both Lasso and Post-Lasso from a cv.plasso object.
#' 
#' @param object \code{\link[plasso]{cv.plasso}} object
#' @param ... Pass generic \code{\link[stats]{coef}} options
#' @param s Determines whether coefficients are extracted for all values of lambda ("all") or only for the optimal lambda ("optimal") according to the specified standard error-rule.
#' @param se_rule Only If equal to zero predictions from CV minimum (default). Negative values go in the direction of smaller
#' models (e.g. se_rule=-1 creates the standard 1SE rule), positive values go in the direction of larger models
#' (e.g. se_rule=1 creates the standard 1SE+ rule). This argument is not used for s="all".
#' @param weights Vector of weights for fitting coefficients.
#' 
#' @return List containing matrices or vector of coefficients for both Lasso and Post-Lasso.
#'
#' @method coef cv.plasso
#'
#' @export 
#' 
coef.cv.plasso = function(object,...,s=c("optimal","all"),se_rule=0,weights=NULL){
  return(predict(object,...,s=s,se_rule=se_rule,type="coefficients",weights=weights))
}


#' Plot of cross-validation curves
#' 
#' @description
#' Plot of cross-validation curves.
#' 
#' @param x \code{\link[plasso]{cv.plasso}} object
#' @param ... Pass generic \code{\link[base]{plot}} options
#' @param legend_pos Legend position. Only considered for joint plot (lass=FALSE).
#' @param legend_size Font size of legend
#' @param lasso If set as True, only the cross-validation curve for the Lasso model is plotted. Default is False.
#' 
#' @method plot cv.plasso
#'
#' @export
#'
plot.cv.plasso = function(x,...,
                          legend_pos=c("bottomright","bottom","bottomleft","left","topleft","top","topright","right","center"),
                          legend_size=0.5,
                          lasso=FALSE) {
  
  plasso = x
  legend_pos = match.arg(legend_pos)
  
  if (!lasso) {
  
    # Standard error of folds
    oneSE_lasso = sqrt(apply(plasso$cv_MSE_lasso,2,var)/plasso$kf)
    oneSE_plasso = sqrt(apply(plasso$cv_MSE_plasso,2,var)/plasso$kf)
    oneSE_lasso[is.na(oneSE_lasso)] = 0
    oneSE_plasso[is.na(oneSE_plasso)] = 0
    
    # Calculate 1SE bands
    lasso_1se_up = plasso$mean_MSE_lasso+oneSE_lasso
    lasso_1se_low = plasso$mean_MSE_lasso-oneSE_lasso
    plasso_1se_up = plasso$mean_MSE_plasso+oneSE_plasso
    plasso_1se_low = plasso$mean_MSE_plasso-oneSE_plasso
    
    # Get ranges for the graph
    xrange = range(log(plasso$lasso_full$lambda))
    yrange = c(max(-1.7e+308,min(lasso_1se_low,plasso_1se_low)),
               min(1.7e+308,max(lasso_1se_up,plasso_1se_up)))
    
    # Plot mean lines
    ylab = "Mean-squared Error"
    graphics::plot(xrange,yrange,type="n",xlab="Log Lambda",ylab=ylab)
    graphics::lines(log(plasso$lasso_full$lambda),plasso$mean_MSE_lasso,lwd=1.5,col="blue")
    graphics::lines(log(plasso$lasso_full$lambda),plasso$mean_MSE_plasso,lwd=1.5,col="red")
    
    # Plot upper and lower 1SE lines
    graphics::lines(log(plasso$lasso_full$lambda),lasso_1se_up,lty=2,lwd=1,col="blue")
    graphics::lines(log(plasso$lasso_full$lambda),lasso_1se_low,lty=2,lwd=1,col="blue")
    graphics::lines(log(plasso$lasso_full$lambda),plasso_1se_up,lty=2,lwd=1,col="red")
    graphics::lines(log(plasso$lasso_full$lambda),plasso_1se_low,lty=2,lwd=1,col="red")
    
    # Show location of minima
    graphics::abline(v=log(plasso$lambda_min_l),lty = 1,col="blue")
    graphics::abline(v=log(plasso$lambda_min_pl),lty = 1,col="red")
    
    # Print legend
    graphics::legend(legend_pos, c("Lasso MSE","Lasso MSE+-1SE","Post-Lasso MSE","Post-Lasso MSE+-1SE","# active coeff."),lty=c(1,2,1,2,1),
                     col=c('blue','blue','red','red','forestgreen'),ncol=1,bty ="n",cex=legend_size)
    
    # Open a new graph for number of coefficients to be written on existing
    graphics::par(new=TRUE)
    graphics::plot(log(plasso$lasso_full$lambda),plasso$lasso_full$df,axes=F,xlab=NA,ylab=NA,cex=1.2,type="l",col="forestgreen",lwd=1.5)
    graphics::axis(side=4)
    graphics::mtext(side=4, line=3, "# active coefficients")
    
  } else if (lasso) {
    
    plot(plasso$lasso_full,...)
    
  }
}


#' Core part for (Post-) Lasso cross-validation
#' 
#' @description
#' \emph{CV_core()} ccontains the core parts of the CV for Lasso and Post-Lasso.
#' 
#' @param x Covariate matrix to be used in CV
#' @param y Vector of outcomes
#' @param w Vector of weight
#' @param cvgroup Categorical with k groups to identify folds
#' @param list List 1:k
#' @param i Number of fold that is used for prediction
#' @param lambda Series of lambdas used
#' @param ... Pass \code{\link[glmnet]{glmnet}} options
#'
#' @return MSE_lasso / MSE_plasso: means squared errors for each lambda.
#'
#' @keywords internal
#'
CV_core = function(x,y,w,cvgroup,list,i,lambda,...) {
  
  # Get estimation and prediction sample for this specific fold
  x_est_cv = x[cvgroup %in% list[-i],]
  y_est_cv = y[cvgroup %in% list[-i]]
  w_est_cv = w[cvgroup %in% list[-i],]
  x_pred_cv = x[cvgroup == list[i],]
  y_pred_cv = y[cvgroup == list[i]]
  w_pred_cv = w[cvgroup == list[i],]
  
  # Normalize the weights to N
  w_est_cv = norm_w_to_n(as.matrix(w_est_cv))
  w_pred_cv = norm_w_to_n(as.matrix(w_pred_cv))
  
  # Estimate Lasso for this fold using the grid of the full sample
  lasso_cv = glmnet::glmnet(x_est_cv, y_est_cv,lambda = lambda,weights=as.vector(w_est_cv),
                    family="gaussian",...)
  coef_lasso_cv = coef(lasso_cv)                                       # Save coefficients at each grid point
  
  # Predicted values with lasso coefficients for prediction sample at each grid
  fit_lasso = predict(lasso_cv,x_pred_cv)
  if (ncol(fit_lasso) != length(lambda)) {
    fit_lasso = cbind(fit_lasso,matrix(NA,nrow(fit_lasso),(length(lambda)-ncol(fit_lasso))))
  }
  
  ## Now calculate predicted values for post Lasso at each grid
  # Initialize first matrix for fitted values
  fit_plasso = matrix(NA,nrow = nrow(fit_lasso),ncol = ncol(fit_lasso))
  
  # Figure out which coefficients are active at each Lasso grid
  act_coef = (coef_lasso_cv != 0)         # boolean matrix
  
  ## Calculate full covariance matrix X'X and X'y once and select only the relevant in the loop below (much faster)
  # First, figure out variables that were active at least once and get only these
  act_once = apply(act_coef,1,any)
  nm_all_act_coef = rownames(act_coef)[act_once]
  if (length(nm_all_act_coef) == 1) {
    warning("No variables selected in one CV, even for lowest lambda, might reconsider choice of penalty term grid")
    fit_plasso = fit_lasso
  } else {
    ## Create the covariate matrix to "manually" calculate the fitted values (much faster than the build in lm command)
    if (sum(abs(coef_lasso_cv[1,])) == 0) {   # Indicates that no intercept was used
      x_all_act = x_est_cv[,nm_all_act_coef]
    }
    else if (sum(abs(coef_lasso_cv[1,])) != 0) {    # Indicates that intercept was used
      x_all_act = add_intercept(x_est_cv[,nm_all_act_coef[2:length(nm_all_act_coef)],drop=F])
      colnames(x_all_act)[1] = "(Intercept)"
      # add intercept also to prediction sample
      x_pred_cv = add_intercept(x_pred_cv[,nm_all_act_coef[2:length(nm_all_act_coef)],drop=F])
      colnames(x_pred_cv)[1] = "(Intercept)"
    }
    else {
      stop("Something strange happens with the intercepts at the Lasso path")
    }
    # Add weights
    x_w = apply(x_all_act,2,`*`,sqrt(w_est_cv))
    y_w = y_est_cv * sqrt(w_est_cv)
    
    # Get X'X
    XtX_all = crossprod(x_w)
    # Get X'y
    Xty_all = crossprod(x_w,y_w)
    
    # Initialize vector to save chosen variable names to figure out whether new variables were added from the last to the next grid
    nm_act_coef_prev = "nobody should call a variable like this"
    
    ## Loop over the grid of Lambdas to get the Post-Lasso predictions
    for (j in 1:length(lambda)) {
      # Get names of active variables at this grid
      nm_act_coef = rownames(act_coef)[act_coef[,j]]
      if (identical(nm_act_coef,character(0)) == TRUE) {
        # print("No variable selected at this grid => no prediction")
        next
      }
      if (identical(nm_act_coef,nm_act_coef_prev) == TRUE & j > 1) {
        # print("Same variables selected as in previous grid => Post-Lasso predictions remain unchanged")
        fit_plasso[,j] = fit_plasso[,(j-1)]
        next
      }
      else {     # Get prediction covariate matrix for that grid
        x_ols_pred = as.matrix(x_pred_cv[,nm_act_coef])
      }
      
      # Get OLS Post-Lasso predictions for this grid point
      fit = fitted_values_cv(XtX_all,Xty_all,x_ols_pred,nm_act_coef)
      if (is.null(fit) & j == 1) {
        fit_plasso[,j] = rep(mean(y),nrow(fit_plasso))
        next
      }
      if (is.null(fit) & j > 1) {
        # cat("\n X'X not invertible at grid",toString(j),": Use last feasible coefficients")
        fit_plasso[,j] = fit_plasso[,(j-1)]
        next
      } else{
        fit_plasso[,j] = fit
      }
      
      # Update current active covariates for the check whether it changes in the next grid
      nm_act_coef_prev = nm_act_coef
    }       # end loop over grids
    
  } # end if at least one var selected
  
  # Matrix with "real" outcome values for each grid
  y_rep = matrix(rep(y_pred_cv,length(lambda)),nrow=nrow(fit_lasso),ncol=ncol(fit_lasso))
  
  # Get RMSE
  SE_lasso = (y_rep - fit_lasso)^2
  SE_plasso = (y_rep - fit_plasso)^2
  if (!is.null(w)) {                                     # Weight errors by "sampling weights"
    SE_lasso = apply(SE_lasso,2,"*",w_pred_cv)
    SE_plasso = apply(SE_plasso,2,"*",w_pred_cv)
  }
  MSE_lasso = apply(SE_lasso,2,mean)
  MSE_plasso = apply(SE_plasso,2,mean)
  
  return(list("MSE_lasso" = MSE_lasso,"MSE_plasso" = MSE_plasso))
}