### This R file contains the functions conducting the simulations of the true 
### PWER presented in section 4 of the paper.

library(ggplot2)  
library(MASS)
library(tidyr)
library(doParallel)
no_cores <- detectCores()-2  
library(foreach)

### Functions that simulate the values of the true PWER in the various different 
### cases of section 4.

# Prevalences: random and based on independent biomarkers, Estimator: MLE, 
# Variance of the response: unknown
f1 <- function(k, m, Nsim, alpha, N, value){
  source("Functions.R")
  p <- runif(m)
  pi.true <- compute_pi(m,p)
  pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
  Sigma <- covmatfct(pi.est, m)
  c.est <- critpwer(m, pi.est, Sigma, alpha, N)
  if(value=="pwer") rval <- pwerfct(m, pi.true, Sigma, c.est, N)
  if(value=="fwer") rval <- max_fwerfct(m, Sigma, c.est, N)
  return(rval)
}

# Prevalences: random and based on correlated biomarkers, Estimator: MLE, 
# Variance of the response: unknown
f2 <- function(k, m, Nsim, alpha, N, value){
  source("Functions.R")
  sigma <- randcorr(m)
  rand_vals <- mvrnorm(n=1,mu=rep(0, m),Sigma=sigma)
  p <- pnorm(rand_vals)
  pi.true <- compute_pi(m,p)
  pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
  Sigma <- covmatfct(pi.est, m)
  c.est <- critpwer(m, pi.est, Sigma, alpha, N)
  if(value=="pwer") rval <- pwerfct(m, pi.true, Sigma, c.est, N)
  if(value=="fwer") rval <- max_fwerfct(m, Sigma, c.est, N)
  return(rval)
}

# Prevalences: equal in all strata, Estimator: MLE, Variance of the response: 
# unknown
f3 <- function(k, m, Nsim, alpha, N, value){
  source("Functions.R")
  pi.true <- matrix(rep(1/(2^m-1),2^m-1), nrow=1, ncol=2^m-1)
  pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
  Sigma <- covmatfct(pi.est, m)
  c.est <- critpwer(m, pi.est, Sigma, alpha, N)
  if(value=="pwer") rval <- pwerfct(m, pi.true, Sigma, c.est, N)
  if(value=="fwer") rval <- max_fwerfct(m, Sigma, c.est, N)
  return(rval)
}

# Prevalences: 0.5 in one stratum and equal in the others, Estimator: MLE, 
# Variance of the response: unknown
f4 <- function(k, m, Nsim, alpha, N, value){
  source("Functions.R")
  pi.true <- matrix(c(0.5, rep(0.5/(2^m-2),2^m-2)), nrow=1, ncol=2^m-1)
  pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
  Sigma <- covmatfct(pi.est, m)
  c.est <- critpwer(m, pi.est, Sigma, alpha, N)
  if(value=="pwer") rval <- pwerfct(m, pi.true, Sigma, c.est, N)
  if(value=="fwer") rval <- max_fwerfct(m, Sigma, c.est, N)
  return(rval)
}

# Prevalences: random and based on independent biomarkers, Estimator: Marginal 
# sum estimator, Variance of the response: unknown
f5 <- function(k, m, Nsim, alpha, N, value){
  source("Functions.R")
  p <- runif(m)
  pi.true <- compute_pi(m,p)
  pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
  pi.est <- marginal_pi(m, pi.est)
  Sigma <- covmatfct(pi.est, m)
  c.est <- critpwer(m, pi.est, Sigma, alpha, N)
  if(value=="pwer") rval <- pwerfct(m, pi.true, Sigma, c.est, N)
  if(value=="fwer") rval <- max_fwerfct(m, Sigma, c.est, N)
  return(rval)
}

# Prevalences: random and based on independent biomarkers, Estimator: MLE, 
# Variance of the response: known
f6 <- function(k, m, Nsim, alpha, N, value){
  source("Functions.R")
  p <- runif(m)
  pi.true <- compute_pi(m,p)
  pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
  Sigma <- covmatfct(pi.est, m)
  c.est <- critpwer_kv(m, pi.est, Sigma, alpha)
  if(value=="pwer") rval <- pwerfct_kv(m, pi.true, Sigma, c.est)
  if(value=="fwer") rval <- max_fwerfct_kv(m, Sigma, c.est)
  return(rval)
}

# Prevalences: random and based on independent biomarkers, Estimator: MLE, 
# Variance of the response: unknown, one single treatment investigated in all
# strata
f7 <- function(k, m, Nsim, alpha, N, value){
  source("Functions.R")
  p <- runif(m)
  pi.true <- compute_pi(m,p)
  pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
  Sigma <- covmatfct_s(pi.est, m)
  c.est <- critpwer(m, pi.est, Sigma, alpha, N)
  if(value=="pwer") rval <- pwerfct(m, pi.true, Sigma, c.est, N)
  if(value=="fwer") rval <- max_fwerfct(m, Sigma, c.est, N)
  return(rval)
}

# Prevalences: random and based on independent biomarkers, Estimator: MLE, 
# Variance of the response: unknown, Treatment allocation numbers in strata: random 
f8 <- function(k, m, Nsim, alpha, N, value){
  source("Functions.R")
  p <- runif(m)
  pi.true <- compute_pi(m,p)
  n <- rmultinom(n=1, size=N, prob=pi.true)
  pi.est <- n/N
  S <- allocation(m,n)
  Sigma <- covmatfct_ra(m, S)
  c.est <- critpwer(m, pi.est, Sigma, alpha, N)
  if(value=="pwer") rval <- pwerfct(m, pi.true, Sigma, c.est, N)
  if(value=="fwer") rval <- max_fwerfct(m, Sigma, c.est, N)
  return(rval)
}

#Calls one of the obove functions Nsim times
sim <- function(m, Nsim, alpha, N, case, value){
  for(i in 1:8) if(case==i) f <- get(paste0("f", i)) 
  cl <- makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  results <- foreach::foreach(seed=1:Nsim) %dopar% {
    set.seed(seed)
    return(f(seed,m,Nsim,alpha,N,value))
  }
  stopCluster(cl)
  return(unlist(results))
}





