#### R codes from the paper: "A shared spatial model for multivariate extreme-valued binary data with non-random missingness", by
#### Xiaoyue Zhao, Lin Zhang, and Dipankar Bandyopadhyay
    
#### The Hamiltonian MC steps for updating $\mu_i$



update_mu_i <- function(sub, mu, accept.mu, L,epsilon, temp.Wbeta, temp.oalpha, full_covariance, inv_cov, xi,y_impute,a_0,b_0,y_m_indicator){
    
    # function to update latent mu with Hamiltonian MCMC
    
    current_q <- q <- matrix(mu[sub,],nrow=1)
    current_p <- p <- matrix(rnorm(length(q)),nrow=1) # independent standard normal variates
    
    
#####  Update with the Leapfrog Method; See Section 2.3 (equations 2.28,29,30) of Radford Neal's "MCMC using Hamiltonian Dynamics", appearing as Chapter 5 in 
#####  Handbook of Markov Chain Monte Carlo. Also, see link: https://arxiv.org/pdf/1206.1901
     
    
    for (k in 1:L){
        
        # Make a half step for momentum at the beginning
        
        p = p - (epsilon[sub]/2) * calculate_gradient(q,sub,temp.Wbeta,temp.oalpha,
                                                      inv_cov,xi,y_impute,a_0,b_0,y_m_indicator)  ##update momentum
        
        # Make a full step for the position
        
        q = q + epsilon[sub] * p
        
        # Make a half step for the momentum, except at end of trajectory
        
        p = p - (epsilon[sub]/2) * calculate_gradient(q,sub,temp.Wbeta,temp.oalpha,
                                                      inv_cov,xi,y_impute,a_0,b_0,y_m_indicator)
        
                 } # close: for (k in 1:L)
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = calculate_U(current_q,sub,temp.Wbeta,temp.oalpha,full_covariance,xi,y_impute,a_0,b_0,y_m_indicator)
    current_K = (current_p%*%t(current_p)) / 2
    proposed_U = calculate_U(q,sub,temp.Wbeta,temp.oalpha,full_covariance,xi,y_impute,a_0,b_0,y_m_indicator)
    proposed_K = (p%*%t(p)) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    # calculate accept.rate for each mu
    
    if (log(runif(1)) < (current_U-proposed_U+current_K-proposed_K)) {
        accept.mu[sub] <- accept.mu[sub] + 1 
        mu[sub,] <- q  # accept
        #cat('move at iteration', r, 'subject', i,'\n')
    }
    
    #return(list(mu.star[r,sub,], accept.mu[sub]))
    return(c(mu[sub,], accept.mu[sub]))
                                                                                                                         } #close update_mu_i function



## Computing the U function
## adjust pgev (0.01, 0.99)

calculate_U <-function(mu_value, sub, temp.Wbeta,temp.oalpha,full_covariance,
                       xi,y_impute,a_0,b_0,y_m_indicator){
    U1 <- -dmvnorm(mu_value-temp.Wbeta[,sub]-temp.oalpha, sigma=full_covariance, log=TRUE)
    
    adjust_gev1 <- sapply(mu_value, function(x) min(max(pgev(-x,shape=xi),0.0001),0.9999))
    U2 <- - y_impute[sub,] %*% log(1-adjust_gev1) - (1 - y_impute[sub,]) %*% log(adjust_gev1)
    
    adjust_pnorm <- sapply(mu_value, function(x) min(max(pnorm(a_0+b_0*x),0.0001),0.9999))
    U3 <- -(y_m_indicator[sub,])%*%log(adjust_pnorm) - (1-y_m_indicator[sub,])%*%log(1 - adjust_pnorm)
    
    U <- U1 + U2 + U3
    return(U)
                                                        }



## grad_U function ##
# gradient = grad_U1 + grad_U2 + grad_U3
# MVN, GEV, Normal



calculate_gradient <- function(mu_value, sub, temp.Wbeta,temp.oalpha,inv_cov,
                               xi,y_impute,a_0,b_0,y_m_indicator){
    
    # grad_U1 <- t(inv_cov%*%mu_value - inv_cov%*%(W[,,sub]%*%beta+omega*alpha))
    grad_U1 <- (mu_value-temp.Wbeta[,sub]-temp.oalpha) %*% inv_cov
    
    adjust_gev <- sapply(mu_value, function(x) min(max(pgev(-x,shape=xi),0.0001),0.9999))
    temp_dgev <- dgev(-mu_value,shape=xi)
    grad_U2 <-  matrix( -y_impute[sub,]*temp_dgev/(1-adjust_gev) + (1- y_impute[sub,])*temp_dgev/adjust_gev,nrow=1)
    
    adjust_pnorm1 <- sapply(mu_value, function(x) min(max(pnorm(a_0+b_0*x),0.0001),0.9999))
    temp_dnorm <- dnorm(a_0+b_0*mu_value)
    grad_U3 <- b_0*matrix(-y_m_indicator[sub,]*temp_dnorm/adjust_pnorm1 + (1-y_m_indicator[sub,])*temp_dnorm/(1-adjust_pnorm1),nrow=1)
    
    grad_U <- grad_U1 + grad_U2 + grad_U3
    return(grad_U)
                                                                 }

