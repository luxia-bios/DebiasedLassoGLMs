###########################################################################
### This R code file provides an example to draw inference for GLMs     ###
### with simulated data. It computes three estimates: the proposed      ###
### refined de-biased lasso estimates, the original de-biased lasso     ###
### estimates, and the MLE.
###                                                                     ###
### Dependencies: glmnet, mvtnorm                                       ###
### Contact: Lu Xia (email: xialu@uw.edu)                               ###
### Date: July 19, 2021                                                 ###
###########################################################################


rm(list=ls())
library(mvtnorm)
library(glmnet)

source("./DBL_GLMs_functions.R")

set.seed(159) # set seed for reproducibility

n_lambda <- 100 # number of lambda values for cross-validation for node-wise lasso
nfold <- 10 # number of folds for cross-validation for node-wise lasso
lambda_ratio <- 0.005
sig_level <- 0.05 # significance level
v <- qnorm(sig_level/2, lower.tail=F)

n <- 500 # number of observations
p <- 100 # number of covariates
struct <- "ar1" # covariance structure for covariates
rho <- 0.7
covmat <- matrix(0, nrow=p, ncol=p) 
if(struct == "indep") {
  covmat <- diag(p)
} else if(struct == "ar1") {
  covmat <- rho^(abs(outer(1:p,1:p,"-")))
} else if(struct == "cs") {
  covmat <- rho*rep(1,p)%*%t(rep(1,p)) + (1-rho)*diag(p)
} else if(struct == "invcs") {
  covmat <- rho*rep(1,p)%*%t(rep(1,p)) + (1-rho)*diag(p)
  covmat <- solve(covmat)
  covmat <- diag(1/sqrt(diag(covmat)))%*%covmat%*%diag(1/sqrt(diag(covmat)))
}
s0 <- 4 # number of additional signals
large_signal <- 1
small_signal <- 0.5
beta_true <- rep(0,p+1) # specify the true betas
beta_true[sample(1:p, size=s0)+1] <- c(rep(small_signal,s0/2), rep(large_signal,s0/2))
beta_true[2] <- 1 # specify the true coefficient for the first real covariate; this varies from 0 to 1.5 in the paper
signal_pos <- which(beta_true != 0) # signal positions in beta_true

### data generation (logistic regression) ###

X_pre <- rmvnorm(n, mean = rep(0,p), sigma = covmat) # without intercept column
X_pre <- ifelse(abs(X_pre)>6, sign(X_pre)*6, X_pre) # truncated covariates
X <- cbind(rep(1,n), X_pre) # added intercept column
y <- rbinom(n, size=1, prob=exp(X%*%beta_true)/(1+exp(X%*%beta_true)))

### cross-validation to obtain lasso estimator ###

cvobj_glmnet <- cv.glmnet(x=X_pre, y=y, family="binomial", alpha=1, standardize=F, intercept=T)
beta_glmnet <- as.vector(coef(glmnet(x=X_pre, y=y, family="binomial", alpha=1, standardize=F,
                                     intercept=T), s=cvobj_glmnet$lambda.min))
beta_glmnet[signal_pos] # check lasso estimates for real signals

### the proposed refined de-biased lasso approach ###

obj_ref_ds <- REF_DS_inf(x=X_pre, y=y, family="binomial", lasso_est=beta_glmnet)
ci_upper_ref_ds <- obj_ref_ds$est + v*obj_ref_ds$se
ci_lower_ref_ds <- obj_ref_ds$est - v*obj_ref_ds$se
cbind(ci_lower_ref_ds, ci_upper_ref_ds) # 95% confidence intervals

### the original de-biased lasso approach ###

obj_orig_ds <- ORIG_DS_inf(x=X_pre, y=y, family="binomial", lasso_est=beta_glmnet)
ci_upper_orig_ds <- obj_orig_ds$est + v*obj_orig_ds$se
ci_lower_orig_ds <- obj_orig_ds$est - v*obj_orig_ds$se
cbind(ci_lower_orig_ds, ci_upper_orig_ds) # 95% confidence intervals

### MLE ###

beta_mle <- coef(glm.fit(x=cbind(rep(1,n),X_pre), y=y, family=binomial()))
mu_mle <- as.vector(1/(1+exp(-X%*%beta_mle)))
neg_ddloglik_mle <- t(X)%*%diag(mu_mle*(1-mu_mle))%*%X/n
se_mle <- sqrt(diag(solve(neg_ddloglik_mle*n)))
pval_mle <- 2*pnorm(abs(beta_mle/se_mle), lower.tail=F)