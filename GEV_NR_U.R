#### R codes from the paper: "A shared spatial model for multivariate extreme-valued binary data with non-random missingness", by
#### Xiaoyue Zhao, Lin Zhang, and Dipankar Bandyopadhyay



## Computing the U function in the Hamiltonian Step## 


calculate_U <-function(mu_value,sub){
U1 <- -dmvnorm(t(mu_value), W[,,sub]%*%beta.star[r,]+omega%*%alpha.star[r,], full_covariance, log=TRUE)
adjust_gev1 <- pgev(-mu_value, scale=1, shape=xi.star[r])
adjust_gev1 <- sapply(1:length(adjust_gev1), function(ii) min(max(adjust_gev1[ii],0.0001),0.9999))

U2 <- -matrix(y_impute[sub,],nrow=1)%*%log(1-adjust_gev1) - matrix((rep(1,n.teeth)-y_impute[sub,]),nrow=1)%*%log(adjust_gev1)

adjust_pnorm <- pnorm(a_0.star[r-1]+b_0.star[r-1]*Z%*%mu_value, 0, 1)
adjust_pnorm <- sapply(1:length(adjust_pnorm), function(jj) min(max(adjust_pnorm[jj],0.0001),0.9999))
U3 <- -(y_m_indicator_teeth[sub,])%*%log(adjust_pnorm) -(1-y_m_indicator_teeth[sub,])%*%log(1 - adjust_pnorm)

U <- U1 + U2 + U3
return(U)
}



## grad_U function ##
# gradient = grad_U1 + grad_U2 + grad_U3
# MVN, GEV, Normal

calculate_gradient <- function(mu_value,sub){
grad_U1 <- t(inv_cov%*%mu_value - inv_cov%*%(W[,,sub]%*%beta.star[r,]+omega%*%alpha.star[r,]))
adjust_gev <- pgev(-mu_value,scale=1, shape=xi.star[r])
adjust_gev <- sapply(1:length(adjust_gev), function(ii) min(max(adjust_gev[ii],0.0001),0.9999))

grad_U2 <- matrix(-y_impute[sub,]*(dgev(-mu_value, scale=1, shape=xi.star[r])/(1-adjust_gev)) + (1- y_impute[sub,])*(dgev(-mu_value, scale=1, shape=xi.star[r])/(adjust_gev)),nrow=1)

adjust_pnorm1 <- pnorm(a_0.star[r-1]+b_0.star[r-1]*Z%*%mu_value, 0, 1)
adjust_pnorm1 <- sapply(1:length(adjust_pnorm1), function(jj) min(max(adjust_pnorm1[jj],0.0001),0.9999))
grad_U3 <- matrix(-y_m_indicator_teeth[sub,]*(dnorm(t(a_0.star[r-1]+b_0.star[r-1]*Z%*%mu_value))/adjust_pnorm1), nrow=1)%*%(b_0.star[r-1]*Z) + matrix((1-y_m_indicator_teeth[sub,])*(dnorm(t(a_0.star[r-1]+b_0.star[r-1]*Z%*%mu_value))/(1-adjust_pnorm1)),nrow=1)%*%(b_0.star[r-1]*Z)
          
grad_U <- grad_U1 + grad_U2 + grad_U3
return(grad_U)
}



