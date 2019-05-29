#### R codes from the paper: "A shared spatial model for multivariate extreme-valued binary data with non-random missingness", by
#### Xiaoyue Zhao, Lin Zhang, and Dipankar Bandyopadhyay



rm(list=ls())

## install some R libraries 
load.lib<-c("evd", "mvtnorm", "msm","truncnorm","coda","doParallel","pscl","plyr","foreach")


install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)


### load the R libraries 
library(MASS)
library(evd)
library(mvtnorm)
library(msm)
library(truncnorm)
library(coda)
library(plyr)
library(foreach)
library(doParallel)
library(pscl)



set.seed(1234)



## Set your favorite working directory
setwd("C:/.../") 



# load("chain1_data_initial.RData")


load("PDData.RData") ## loading the dataset

source("GEV_NR_U.R") ## R codes for the Hamiltonian MC steps



CORE <- 4  ### Here, we are using 4 cores to run the program 


nn <- 10000
initial_beta <- beta.save[nn,]
initial_alpha <- alpha.save[nn]
initial_sigma2 <- sigma2.save[nn]
initial_mu <- mu.save[nn,,]
initial_rho <- rho.save[nn]
initial_xi <- xi.save[nn]
initial_a_0 <- a_0.save[nn]
initial_b_0 <- b_0.save[nn]
epsilon <- epsilon_hist[,50]
 



NB <- 20000
N <-  10000


# epsilon <- c(rep(0.05, n.subs)) # c(rep(0.005, n.subs))

stepsize <- 0.005 # 0.0005
L <- 150 

n_check_accept <- 100
n_check_trace <- 1000



### Arrays to save parameters

beta.save <- array(NA, c(N, p1)) 
alpha.save <- array(NA, c(N, p2)) 
a_0.save <- rep(NA, N)
b_0.save <- rep(NA, N)
sigma2.save <- rep(NA, N)
rho.save <- rep(NA, N)
xi.save <- rep(NA, N) 
mu.save <- array(NA, c(N, n.subs, n.teeth))
# y_latent.save <- array(NA, c(N, n.subs, n.teeth)) 



### Initial values

initial_a_0 <- -1 
initial_xi <- -2 
initial_rho <- 0.975 
initial_sigma2 <- 0.25

beta <- initial_beta
alpha <- initial_alpha
sigma2 <- initial_sigma2
mu <- initial_mu
rho <- initial_rho
xi <- initial_xi
a_0 <- initial_a_0
b_0 <-initial_b_0
y_latent <- matrix(NA,n.subs,n.teeth) # initial_y_i0star

temp.Wbeta <- sapply(alply(W,3), function(x1) x1 %*% beta, simplify='cbind')
temp.oalpha <- omega*alpha
V_rho = diag(m) - rho*ADJ

y_m_indicator <- is.na(CAL2) # missingness indicator, 1==missing, 0==present



### Hyperparameters
a_sigma2 <- b_sigma2 <- 0.1
# a_sigma2 <- b_sigma2 <- n.subs*n.teeth*5

alpha_variance_tun <- diag(100,p2,p2)
inv_alpha_variance_tun <- solve(alpha_variance_tun)
beta_variance_tun <-  diag(100,p1,p1)
inv_beta_variance_tun <- solve(beta_variance_tun)
beta_a <- rep(0, p1)
alpha_a <- rep(0, p2)

xi.step <- 0.1
xi_bounds <- c(-10,10)
rho_a <- 1
rho_b <- 0.3 # concentrate at 1



### Record acceptance rates for MH steps
accept.xi <- 0
accept.rho <- 0
accept.mu <- rep(0,n.subs)
accept_rate_mu <- rep(NA, n.subs)
accept_rate_hist <- epsilon_hist <- integer(0)


##### MCMC interations, and parameter updatings


# cl= makeCluster(CORE)

registerDoParallel(CORE)

for (r in 1:(NB+N)){
    
    #impute missingness
    y_impute <- CAL2
    
    temp_p <- sapply(mu[y_m_indicator], function(x) 1-pgev(-x,shape=xi))
    y_impute[y_m_indicator] <- sapply(temp_p, function(x) rbinom(1,size=1,prob=x))
    
    
## update sigma2 #debugged
    if (TRUE){
        # if (FALSE) {
        
        # temp.Wbeta <- sapply(alply(W,3), function(x1) x1 %*% beta + omega*alpha, simplify='cbind' )
        temp.u <- mu-t(temp.Wbeta)-matrix(temp.oalpha,n.subs,n.teeth,byrow=T)
        
        temp.sigma2 <- sum(apply(temp.u,1,function(x) t(x)%*%V_rho %*% x) )/2
        sigma2 <- 1/rgamma(1, (n.subs*n.teeth)/2 + a_sigma2, b_sigma2 + temp.sigma2)
        
    } 
    
    
## update_beta #debugged
    if (TRUE){
        # if (FALSE){
        
        temp.G <- sapply(alply(W,3),function(x) t(x)%*%V_rho%*%x,simplify='array')
        temp.G <- apply(temp.G,1:2,sum)/sigma2
        
        temp1 <- mapply(function(x1,x2) (x1-temp.oalpha)%*%V_rho%*%x2, alply(mu,1),alply(W,3))
        temp1 <- rowSums(temp1)/sigma2
        
        G <- solve(inv_beta_variance_tun + temp.G)
        g <- c(beta_a%*%inv_beta_variance_tun + temp1)
        
        beta <- mvrnorm(n=1,mu=G%*%g, Sigma=G) # mean, Sigma
        temp.Wbeta <- sapply(alply(W,3), function(x1) x1 %*% beta, simplify='cbind' )
    } 
    
    
    
## update alpha #debugged
    if (TRUE){
        # if (FALSE){
        
        temp_H <- t(omega)%*%V_rho%*%omega/sigma2
        
        # temp.Wbeta <- mapply(function(x1,x2) x1 %*% x2, alply(W,3), rep(list(beta),n.subs) )
        temp2 <- sum( (mu-t(temp.Wbeta)) %*% V_rho %*% omega /sigma2 )
        
        H <- solve(inv_alpha_variance_tun + n.subs*(temp_H))
        h <- alpha_a %*% inv_alpha_variance_tun + temp2 
        
        alpha <- mvrnorm(1, H%*%h, H ) #mean, Sigma
        temp.oalpha <- omega*alpha
    } 
    
    
    
## update_rho #low acceptance rate
    if (TRUE){
        # if (FALSE){
        
        cons <- 1000
        
        rho.can <- min(0.9999, max(0.95, rbeta(1, cons*rho, cons*(1-rho))))
        V_rho.can <- diag(m) - rho.can*ADJ
        
        temp_cov <- chol2inv(chol(V_rho.can))*sigma2
        temp5 <- dmvnorm(mu-t(temp.Wbeta),temp.oalpha,temp_cov,log=TRUE)
        
        temp_cov_cur <- chol2inv(chol(V_rho))*sigma2
        temp6<- dmvnorm(mu-t(temp.Wbeta),temp.oalpha,temp_cov_cur,log=TRUE)
        
        rho.frac.numerator <- dbeta(rho.can, rho_a, rho_b, log=TRUE) + sum(temp5) + dbeta(rho,cons*rho.can, cons*(1-rho.can),log=TRUE)  #distribution at candidate, proposal at (current|candidate)
        rho.frac.denominator <- dbeta(rho, rho_a, rho_b, log=TRUE) + sum(temp6) + dbeta(rho.can,cons*rho, cons*(1-rho),log=TRUE) #distribution at current r-1,proposal at (candidate|current)
        rho.ratio <- rho.frac.numerator - rho.frac.denominator
        
        if( log(runif(1)) < rho.ratio){
            accept.rho <- accept.rho+1
            rho <- rho.can
            V_rho = diag(m) - rho*ADJ
        } 
    } 
    
    
 ## update xi #debugged
    if (TRUE){
        # if (FALSE){
        xi.ss <- xi.step
        xi.can <- min(max( (runif(1, min=-xi.ss, max=xi.ss)+ xi) , xi_bounds[1]), xi_bounds[2])
        
        ind1 <- which(y_impute==1)
        ind0 <- which(y_impute==0)
        
        xi.frac.numerator <- c(sapply(mu[ind1], function(x) max(-800,log(pgev(-x,shape=xi.can,lower.tail=FALSE)))), 
                               sapply(mu[ind0], function(x) max(-800,log(pgev(-x,shape=xi.can,lower.tail=TRUE)))) )
        
        xi.frac.denominator <- c(sapply(mu[ind1], function(x) max(-800,log(pgev(-x,shape=xi,lower.tail=FALSE)))),
                                 sapply(mu[ind0], function(x) max(-800,log(pgev(-x,shape=xi,lower.tail=TRUE)))) )
        
        xi.ratio <- sum(xi.frac.numerator) - sum(xi.frac.denominator)
        
        if(log(runif(1)) < xi.ratio){
            accept.xi <- accept.xi+1
            xi <- xi.can
        } 
    } 
    
    
 ## update mu
    if (TRUE){
        # if (FALSE){
        
        mu.old <- mu
        
        full_covariance <- chol2inv(chol(V_rho))*sigma2 ##covariance matrix for /mu, n.subs*n.teeth
        inv_cov <- V_rho/sigma2
        
        compute <- foreach(i=1:n.subs, .combine='rbind', .packages=c("evd", "mvtnorm")) %dopar% {
            update_mu_i(i,mu.old,accept.mu,L,epsilon, temp.Wbeta, temp.oalpha, 
                        full_covariance,inv_cov,xi,y_impute,a_0,b_0,y_m_indicator)
        }
        mu <- matrix(compute[,1:n.teeth],nrow=n.subs)
        accept.mu <- compute[,n.teeth+1]
        
        # print(dim(mu))
        
        ## check mu acceptance rate, 
        #if stepsize is too small, the acceptance rate will be too high
        
        if (r <= NB & r %% n_check_accept==0){
            # if (r %% n_check_accept==0){
            
            accept_rate_mu <- accept.mu/n_check_accept
            
            cat('acceptance rate for mu at iteration=', r, '\n')
            print(data.frame(epsilon=epsilon,accept_rate=accept_rate_mu))
            
            epsilon.new <- epsilon + (accept_rate_mu > 0.8)*stepsize - (accept_rate_mu < 0.5)*stepsize
            epsilon <- sapply(1:n.subs, function(ii) if(epsilon.new[ii] < stepsize) epsilon[ii]/2 else epsilon.new[ii])
            
            epsilon_hist <- cbind(epsilon_hist,epsilon)
            accept_rate_hist <- cbind(accept_rate_hist, accept_rate_mu)
            
            accept.mu <- rep(0,n.subs)
            
            
        } #close if (r %% n_check_accept==0)
    } 
    
    
## update y_latent, a0, b0 (for missingnes) #debugged 
    
    if (TRUE){
        # if (FALSE){
        
        ind.m <- which(y_m_indicator==1)  # missing teeth
        ind.nm <- which(y_m_indicator==0)  # not missing teeth
        
        y_latent[ind.m] <- sapply(mu[ind.m], function(x) rtnorm(1,lower=0,mean=a_0+b_0*x))
        y_latent[ind.nm] <- sapply(mu[ind.nm], function(x) rtnorm(1,upper=0,mean=a_0+b_0*x))
        
        
        temp_1mu <- cbind(1,c(mu))
        bigm <- t(temp_1mu) %*% temp_1mu
        
        smallm <- c(sum(y_latent), sum(y_latent*mu))  # matrix(y_latent,nrow=1) %*% temp_1mu
        
        ab_Sigma <- solve(bigm)
        ab_mean <- ab_Sigma %*% smallm
        
        parm <- rmvnorm(1, ab_mean, ab_Sigma)
        a_0 <- parm[1]
        b_0 <- parm[2]
    } 
    
    
    if(r > NB) {
        rr = r-NB
        beta.save[rr,] <- beta
        alpha.save[rr,] <- alpha
        a_0.save[rr] <- a_0
        b_0.save[rr] <- b_0
        sigma2.save[rr] <- sigma2
        rho.save[rr] <- rho
        xi.save[rr] <- xi
        mu.save[rr,,] <- mu
    }
    
    # if (r > NB & r %% n_check_trace==0) 
    #     save.image(file = "chain1_GEV_nonRandom190204.rda")


    if (r %% 100 == 0) print(r)
    
} #end of mcmc iterations

# stopCluster(cl)



print(data.frame(epsilon=epsilon,accept_rate=accept_rate_mu))

accept_xi_rho=c(accept.xi,accept.rho)/(NB+N)




# res <- list(beta=beta.save,alpha=alpha.save,a0=a_0.save,b0=b_0.save,sigma2=sigma2.save,
#             rho=rho.save,xi=xi.save,mu=mu.save,epsilon=epsilon_hist,accept_mu=accept_rate_hist,
#             accept_xi_rho=accept_xi_rho)


save(beta.save,alpha.save,sigma2.save,rho.save,xi.save,mu.save,a_0.save,b_0.save, epsilon_hist,accept_rate_hist,accept_xi_rho, file="yourfilename.rda")











