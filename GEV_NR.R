#### R codes from the paper: "A shared spatial model for multivariate extreme-valued binary data with non-random missingness", by
#### Xiaoyue Zhao, Lin Zhang, and Dipankar Bandyopadhyay



rm(list=ls())

## install some R libraries 
load.lib<-c("evd", "mvtnorm", "msm","truncnorm","coda","doParallel","pscl")


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
library(doParallel)
library(pscl)


## Set your own seed 

set.seed(1234)


## Set your favorite working directory
setwd("C:/.../") 

load("PDData.RData")
source("GEV_NR_U.R")
CORE <- 20 ## In this program, we are using a multicore server to run the program 

NB <- 20000
N <-  10000


temp_sub <-seq(1, n.subs)

epsilon <- c(rep(0.002, n.subs))
stepsize <- 0.0001
L <- 150 

n_check_accept <- 600
n_check_trace <- 600


#initialize parameters

G <- array(NA, c(NB+N,p1, p1))
g <- array(NA, c(NB+N,p1, 1))
H <- array(NA, c(NB+N,p2, p2))
h <- array(NA, c(NB+N,p2, 1))
beta.star <- array(NA, c(NB+N, p1)) 
alpha.star <- array(NA, c(NB+N, p2)) 
a_0.star <- rep(NA, NB+N)
b_0.star <- rep(NA, NB+N)
sigma2.star <- rep(NA, NB+N)
rho.star <- rep(NA, NB+N)
xi.star <- rep(NA, NB+N) #the shape parameter 
mu.star <- array(NA, c(NB+N, n.subs, n.teeth))
y_latent.star <- array(NA, c(NB+N, n.subs, n.teeth)) 



#initial values

initial_a_0 <- -1 #########
initial_xi <- -2  #########
initial_rho <- 0.975 ########
beta.star[1,] <- initial_beta
alpha.star[1,] <- initial_alpha
sigma2.star[1] <- initial_sigma2
mu.star[1, , ] <- initial_mu
rho.star[1] <- initial_rho
xi.star[1] <- initial_xi
a_0.star[1] <- initial_a_0
b_0.star[1] <-initial_b_0
y_latent.star[1, ,] <- initial_y_i0star


#hyperparameters
a_sigma2 <- 0.1
b_sigma2 <- 0.1


#tuning parameters for the proposal (priors)

alpha_variance_tun <- diag(100,p2,p2)
inv_alpha_variance_tun <- solve(alpha_variance_tun)
beta_variance_tun <-  diag(100,p1,p1)
inv_beta_variance_tun <- solve(beta_variance_tun)
beta_a <-c(rep(0, p1))
alpha_a <-c(rep(0, p2))
xi_variance_tun <- 0.1
mu_prop.sd = 0.1
rho_a <- 1
rho_b <- 0.3 # concentrate at 1
accept.xi <- 0
accept.rho <- 0
accept.mu <- c(rep(0,n.subs))
accept_rate_mu <- c(rep(NA, n.subs))

y_m_indicator <- t(sapply(1:n.subs, function(i) sapply(1:n.teeth, function(s) ifelse(is.na(CAL2[i,s]),1,0)))) #missingness indicator, 1==missing, 0==present
y_m_indicator_teeth <- y_m_indicator #missingness at teeth level

temp7Sum <-c()
temp8Sum <-c()
temp9Sum <- c()
temp10Sum <- c()


#Indicator function

Indicator<-function(x)
    { 
        ifelse(x==1,1,0)
    }


accept_rate_hist <- matrix(NA, ncol=n.subs, nrow=(N+NB)/n_check_accept)
epsilon_hist <- matrix(NA, ncol=n.subs, nrow=(N+NB)/n_check_accept)

update_mu_i <- function(sub, mu.star,accept.mu, r){

current_q <- t(mu.star[r-1,sub,])
  q <- t(as.matrix(mu.star[r-1, sub,]) )
  p <- t(as.matrix(c(rnorm(n.teeth,0,1) ))) # independent standard normal variates
  current_p <- p




 #####  Update with the Leapfrog Method. 
 #####  See Section 2.3 of Radford Neal's "MCMC using Hamiltonian Dynamics", appearing as Chapter 5 in Handbook of Markov Chain Monte Carlo

 for (k in 1:L){
   
  # Make a half-step for momentum at the beginning

      p = p - (epsilon[sub]/2) * calculate_gradient(mu_value=t(q),sub)  ##update momentum
 
  # Make a full step for the position

      q = q + epsilon[sub] * p

  # Make a half step for the momentum, except at end of trajectory

      p = p - (epsilon[sub]/2) * calculate_gradient(mu_value=t(q),sub)
    
              } # close: for (k in 1:L)
  
  # Negate momentum at end of trajectory to make the proposal symmetric

     p = -p

  # Evaluate potential and kinetic energies at start and end of trajectory

     current_U = calculate_U(t(current_q),sub)
     current_K = (current_p%*%t(current_p)) / 2
     proposed_U = calculate_U(t(q),sub)
     proposed_K = (p%*%t(p)) / 2

  # Accept or Reject the state at end of trajectory, returning either
  # the position at the end of the trajectory, or the initial position

  # calculate accept.rate for each mu

  if (log(runif(1)) < (current_U-proposed_U+current_K-proposed_K))
  {
    accept.mu[sub] <- accept.mu[sub] + 1 
    mu.star[r,sub,] <- q  # accept
    #cat('move at iteration', r, 'subject', i,'\n')
  }
  
  else
  
  {
    mu.star[r,sub,]<- current_q  #reject, mu.star[r,i,] = mu.star[r-1,i,]
    #cat('not moving at iteration', r, 'subject', i,'\n')
  }

#return(list(mu.star[r,sub,], accept.mu[sub]))
return(c(mu.star[r,sub,], accept.mu[sub]))
} #close update_mu_i function


### Functions to update $\xi$ 

if (TRUE){

part1 <- function(i,teeth){ #if not missing, calculate; if missing, =0
temp7 <- sapply(teeth, function(ss) ifelse(y_m_indicator[i,ss]==0, log((1-min(max(pgev(-mu.star[r-1,i,ss], scale=1, shape=xi.can),0.01), 0.99 ))^(CAL2[i,ss])),0))
 return(sum(unlist(temp7)))
} 


part2 <- function(i,teeth){ 
temp8 <- sapply(teeth, function(ss) ifelse(y_m_indicator[i,ss]==0, log(min(max(pgev(-mu.star[r-1,i,ss], scale=1, shape=xi.can),0.01), 0.99)^(1-CAL2[i,ss])),0))
return(sum(unlist(temp8)))
}

part3 <- function(i,teeth){ 
temp9 <- sapply(teeth, function(ss) ifelse(y_m_indicator[i,ss]==0, log((1-min(max(pgev(-mu.star[r-1,i,ss], scale=1, shape=xi.star[r-1]), 0.01), 0.99 ))^(CAL2[i,ss])), 0))
return(sum(unlist(temp9)))
}

part4 <- function(i,teeth){ 
temp10 <- sapply(teeth, function(ss) ifelse(y_m_indicator[i,ss]==0, log(min(max(pgev(-mu.star[r-1,i,ss], scale=1, shape=xi.star[r-1]), 0.01), 0.99)^(1-CAL2[i,ss])), 0))
return(sum(unlist(temp10)))
} 

} #close if (TRUE)

#ptm<- proc.time()

#cl= makeCluster(CORE)
registerDoParallel(CORE)
#getDoParWorkers()

# Updating
for (r in 2:(NB+N)){

#### Impute missingness
 y_impute <- CAL2
 for (i in 1:n.subs){ 
   for (ss in 1:n.teeth){
   if (is.na(CAL2[i, ss])){
 temp_p <- 1-pgev(-mu.star[r-1,i,ss],loc=0, scale=1, shape=xi.star[r-1],lower.tail=TRUE)
 y_impute[i,ss]<- rbinom(1, size=1, prob=temp_p)}
   }
  }

### Update $\sigma^2$ #debugged

if (TRUE){
#if (FALSE) { 
temp.sigma2 <- 0
for (i in 1:n.subs){
 daed <- 0.5*(t(mu.star[r-1,i,]-(W[,,i]%*%beta.star[r-1,]+ omega%*%alpha.star[r-1,])))%*%(diag(m)-rho.star[r-1]*ADJ)%*%(mu.star[r-1,i,]-(W[, , i]%*%beta.star[r-1,]+ omega%*%alpha.star[r-1,]))
 temp.sigma2 <- temp.sigma2 + daed 
}
sigma2.star[r] <- 1/rgamma(1, (n.subs*n.teeth)/2 + a_sigma2, b_sigma2 + temp.sigma2)
 } else
sigma2.star[r] <- initial_sigma2 #


### Update $\beta$  #debugged

if (TRUE){
#if (FALSE){
temp.G <- matrix(0,nrow=p1,ncol=p1)
temp1 <- c(rep(0,p1))
sumx <- matrix(0, nrow=n.teeth, ncol=p1)
 for (i in 1:n.subs){
 ewew <- t(W[, ,i])%*%((diag(m)-rho.star[r-1]*ADJ)*(1/sigma2.star[r]))%*%W[,,i]
 temp.G <- temp.G + ewew
 aapp <- mu.star[r-1,i,]%*%((diag(m)-rho.star[r-1]*ADJ)*(1/sigma2.star[r]))%*%W[,,i]
 temp1 <- temp1 + aapp
 zzx <- W[,,i]
 sumx <- sumx + zzx 
}

G[r-1,,] <- solve(inv_beta_variance_tun + temp.G)
temp2 <- alpha.star[r-1,]%*%t(omega)%*%((diag(m)-rho.star[r-1]*ADJ)*(1/sigma2.star[r]))%*%sumx
g[r-1,,] <- t( c(t(beta_a)%*%(inv_beta_variance_tun) + temp1- temp2))
beta.star[r,] <- mvrnorm(n=1,mu=G[r-1,,]%*%g[r-1,,], Sigma=G[r-1,,]) # mean, Sigma
} else
 beta.star[r,] <- initial_beta


### Update $\alpha$ #debugged

if (TRUE){
#if (FALSE){
temp_H <- t(omega)%*%((diag(m)-rho.star[r-1]*ADJ)*(1/sigma2.star[r]))%*%omega
H[r-1,,] <- solve(inv_alpha_variance_tun + n.subs*(temp_H))

temp3 <- c(rep(0,p2))
temp4 <- c(rep(0,p2))
for (i in 1:n.subs){
 wesd <-t(mu.star[r-1,i,])%*% ((diag(m)-rho.star[r-1]*ADJ)*(1/sigma2.star[r]))%*%omega
 temp3 <- temp3 + wesd
 zzx<- t(beta.star[r,])%*%t(W[,,i])%*%((diag(m)-rho.star[r-1]*ADJ)*(1/sigma2.star[r]))%*%omega
 temp4 <- temp4 + zzx
}
h[r-1,,] <- t( t(alpha_a)%*%(inv_alpha_variance_tun) + temp3- temp4)
alpha.star[r,] <- mvrnorm(1, H[r-1,,]%*%h[r-1,,], H[r-1,,] ) #mean, Sigma
} else
alpha.star[r,] <- initial_alpha


### Update $\rho$  #low acceptance rate

if (TRUE){
#if (FALSE){
rho.can <- min(0.9999, max(0.95, rbeta(1, 50*rho.star[r-1], 50*(1-rho.star[r-1]))))
temp_cov <- chol2inv(chol(diag(m)-rho.can*ADJ))*sigma2.star[r]
temp5<- sapply(1:n.subs, function(i) dmvnorm(mu.star[r-1,i,], W[,,i]%*%beta.star[r,]+omega%*%alpha.star[r,], temp_cov, log=TRUE))
rho.frac.numerator <- dbeta(rho.can, rho_a, rho_b, log=TRUE) + sum(temp5) + dbeta(rho.star[r-1],50*rho.can, 50*(1-rho.can),log=TRUE)  #distribution at candidate, proposal at (current|candidate)

temp_cov_cur <- chol2inv(chol(diag(m)-rho.star[r-1]*ADJ))*sigma2.star[r]
temp6<- sapply(1:n.subs, function(i) dmvnorm(mu.star[r-1,i,], W[,,i]%*%beta.star[r,]+omega%*%alpha.star[r,], temp_cov_cur, log=TRUE))
rho.frac.denominator <- dbeta(rho.star[r-1], rho_a, rho_b, log=TRUE) + sum(temp6) + dbeta(rho.can,50*rho.star[r-1], 50*(1-rho.star[r-1]),log=TRUE) #distribution at current r-1,proposal at (candidate|current)
rho.fraction <- rho.frac.numerator - rho.frac.denominator
rho.MH.fraction <- min (log(1), rho.fraction)
if( log(runif(1))<rho.MH.fraction){
accept.rho <- accept.rho+1
rho.star[r] <- rho.can
} else{
rho.star[r]<- rho.star[r-1]
}
} else
rho.star[r] <- initial_rho


### Update $\xi$  #debugged

 if (TRUE){
#if (FALSE){
 xi.can <- min(max( (runif(1, min=-0.1, max=0.1)+ xi.star[r-1]) , -3), 3)

xi.frac.numerator <- dunif(xi.can, min=-3, max=3, log=TRUE)
for (i in 1:n.subs){  # sum over 1:n.subs  
ind1 = which(CAL2[i,]==1)
ind0 = which(CAL2[i,]==0)
xi.frac.numerator <- xi.frac.numerator + part1(i,ind1) + part2(i,ind0)
 } 

xi.frac.denominator <- dunif(xi.star[r-1], min=-3, max=3, log=TRUE)
for (i in 1:n.subs){
ind1 = which(CAL2[i,]==1)
ind0 = which(CAL2[i,]==0)
xi.frac.denominator <- xi.frac.denominator + part3(i,ind1) + part4(i,ind0)
 }

xi.fraction <- xi.frac.numerator - xi.frac.denominator
xi.MH.fraction <- min (log(1), xi.fraction)
if(log(runif(1)) < xi.MH.fraction){
accept.xi <- accept.xi+1
xi.star[r] <- xi.can
 } else 
xi.star[r] <- xi.star[r-1]
# end of if(true/false) statement 
} else
xi.star[r]<- initial_xi


### Update $\mu$ 

 if (TRUE){
#if (FALSE){

full_covariance <- chol2inv(chol(diag(m)-rho.star[r]*ADJ))*sigma2.star[r] ##covariance matrix for /mu, n.subs*n.teeth
inv_cov <- chol2inv(chol(full_covariance))
compute <- foreach(i=1:length(temp_sub), .combine=cbind, .packages=c("evd", "mvtnorm")) %dopar% update_mu_i(sub=temp_sub[i], mu.star=mu.star, accept.mu=accept.mu, r=r)
mu.star[r,,] <- t(compute[1:n.teeth,])
accept.mu <- compute[n.teeth+1,]


### Check acceptance rate of $\mu$;  
#if stepsize is too small, the acceptance rate will be too high

 if (r %% n_check_accept==0){
    for ( i in 1:n.subs){
 accept_rate_mu[i] <- accept.mu[i]/n_check_accept
  
 if (accept_rate_mu[i]>0.98){
   epsilon[i] <- epsilon[i] + stepsize
  }
  if (accept_rate_mu[i] < 0.75){ #get even smaller steps
   epsilon[i] <- epsilon[i] - stepsize
  } 
 epsilon_hist[(r/n_check_accept),i]<- epsilon[i]
 accept_rate_hist[(r/n_check_accept),i] <- accept_rate_mu[i]
  } # close for (1 in 1:n.subs) 
   accept.mu <- c(rep(0,n.subs))

cat('acceptance rate for mu at iteration=', r, '\n')
print(accept_rate_hist[(r/n_check_accept),])
print(epsilon_hist[(r/n_check_accept),])

    } #close if (r %% n_check_accept==0)

} else
mu.star[r, ,] <- initial_mu


### Update in parallel (latent) $y$ #debugged 

if (TRUE){
#if (FALSE){
    for (i in 1:n.subs){
       for (t in 1:n.teeth){
    if(y_m_indicator_teeth[i,t]==1)
        y_latent.star[r,i,t] <- rtnorm(1,lower=0, upper=Inf, mean= a_0.star[r-1] + t(Z[t,])%*%mu.star[r,i,]*b_0.star[r-1], sd=1) else
        y_latent.star[r,i,t] <- rtnorm(1, lower= -Inf, upper=0, mean=a_0.star[r-1] + t(Z[t,])%*%mu.star[r,i,]*b_0.star[r-1], sd=1)
      } }
    } else
  y_latent.star[r,,] <- initial_y_i0star #50*7
 
 

### Update $a_0$, $b_0$  #### debugged
 
if (TRUE){
#if (FALSE){
 bigm <- matrix(0,2,2)
 smallm<- c(0,0)
 
##calculate the design matrix for n subjects
 
 Dmatrix <- array(NA, c(n.subs,n.teeth,2))
 for (i in 1:n.subs){
 Dmatrix[i,,] <- t(sapply(1:n.teeth, function(sfr) c(1, t(Z[sfr, ])%*%mu.star[r,i,])))
 }

for(i in 1:n.subs) {
  for (t in 1:n.teeth){
     bigm <- bigm +  Dmatrix[i,t,]%*%t(Dmatrix[i,t,]) 
     smallm <- smallm + y_latent.star[r,i,t]*t(Dmatrix[i,t,])  
      }
    }

ab_mean <- solve(bigm)%*%t(smallm)
ab_Sigma <- solve(bigm)
parm <- rmvnorm(1, ab_mean, ab_Sigma)
a_0.star[r] <- parm[1]
b_0.star[r] <- parm[2]

} else {
a_0.star[r] <- initial_a_0
b_0.star[r] <- initial_b_0
}

if (r %% n_check_trace==0){
save(list = ls(all=TRUE), file = "Choose_a_name.RData")
 }
} #iteration close

