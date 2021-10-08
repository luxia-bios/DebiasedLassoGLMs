

###########################################################################
### This R code file provides two functions to draw inference for GLMs. ###
### The function "REF_DS_inf" implements the proposed refined de-biased ###
### lasso approach by directly inverting the Hessian matrix. The other  ###
### function "ORIG_DS_inf" implements the original de-biased lasso by   ###
### van de Geer et al. (2014).                                          ###
###                                                                     ###
### Dependencies: glmnet                                                ###
### Contact: Lu Xia (email: xialu@uw.edu)                               ###
### Date: July 19, 2021                                                 ###
###########################################################################

REF_DS_inf <- function(x, y, family, lasso_est) {
  nn <- length(y)
  X <- cbind(rep(1, nrow(x)), x)
  if(ncol(X) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "binomial") {
    mu <- as.vector(exp(X%*%lasso_est)/(1+exp(X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu*(1-mu))%*%X/nn
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
  } else {
    stop("Input family is not supported.")
  }
  
  theta_inv <- solve(neg_ddloglik_glmnet)
  b_hat_inv <- as.vector(lasso_est - theta_inv%*%neg_dloglik_glmnet)
  se_inv <- sqrt(diag(theta_inv))/sqrt(nn)
  pval_inv <- 2*pnorm(abs(b_hat_inv/se_inv), lower.tail=F)
  
  return(list(est=b_hat_inv, se=se_inv, pvalue=pval_inv, theta=theta_inv))
}


ORIG_DS_inf <- function(x, y, family, lasso_est, nfold=5, n_lambda=100, lambda_ratio=0.005) {
  nn <- length(y)
  pp <- ncol(x)
  X <- cbind(rep(1, nrow(x)), x)
  if(ncol(X) != length(lasso_est)) {
    stop("The length of lasso_est is incompatible with the covariate matrix.")
  }
  if(family == "binomial") {
    mu <- as.vector(exp(X%*%lasso_est)/(1+exp(X%*%lasso_est)))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu*(1-mu))%*%X/nn
    C_glmnet <- sqrt(diag(mu*(1-mu))/nn)%*%X
  } else if(family == "poisson") {
    mu <- as.vector(exp(X%*%lasso_est))
    neg_dloglik_glmnet <- 0 - as.vector(t(X)%*%(y-mu))/nn
    neg_ddloglik_glmnet <- t(X)%*%diag(mu)%*%X/nn
    C_glmnet <- sqrt(diag(mu)/nn)%*%X
  } else {
    stop("Input family is not supported.")
  }
  
  theta_glmnet <- diag(pp+1)
  tau_glmnet <- rep(NA, pp+1)
  for(j in 1:(pp+1)) { # for: nodewise lasso
    current_x <- sqrt(nn)*C_glmnet[,-j]
    current_y <- sqrt(nn)*as.vector(C_glmnet[,j])
    lam_max <- max(abs(t(current_x)%*%current_y)/nn)
    lam_min <- lam_max*lambda_ratio
    lam_seq <- exp(seq(from=log(lam_max), to=log(lam_min), length.out=n_lambda))
    gamma_j_glmnet <- cv.glmnet(x=current_x, y=current_y,
                                family="gaussian", alpha=1, standardize=F, intercept=F,
                                nfolds=nfold, lambda=lam_seq)
    gamma_j_glmnet <- as.vector(glmnet(x=sqrt(n)*C_glmnet[,-j], y=sqrt(n)*as.vector(C_glmnet[,j]),
                                       family="gaussian", alpha=1, standardize=F, intercept=F,
                                       lambda=gamma_j_glmnet$lambda.min)$beta)
    theta_glmnet[j,-j] <- (-1)*t(gamma_j_glmnet)
    tau_glmnet[j] <- as.numeric(neg_ddloglik_glmnet[j,j]-neg_ddloglik_glmnet[j,-j]%*%gamma_j_glmnet)
  } # end for: nodewise lasso
  theta_glmnet <- diag(1/tau_glmnet)%*%theta_glmnet 
  
  b_hat_nw <- as.vector(beta_glmnet - theta_glmnet%*%neg_dloglik_glmnet)
  se_nw <- sqrt(diag(theta_glmnet%*%neg_ddloglik_glmnet%*%t(theta_glmnet)))/sqrt(nn)
  pval_nw <- 2*pnorm(abs(b_hat_nw/se_nw), lower.tail=F)
  
  return(list(est=b_hat_nw, se=se_nw, pvalue=pval_nw, theta=theta_glmnet))
}
