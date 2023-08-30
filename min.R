library(ggplot2)  
library(MASS)
library(tidyr)
library(doParallel)
no_cores <- detectCores()-2
library(gridExtra)

sim_min <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    p <- runif(m,0,0.2)
    pi.true <- compute_pi(m,p)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    
    # minimal prevalence
    min <- 1/(2^(m+1)-2)
    
    # test if one prevalence is 0
    counter = 0
    for(i in 1:length(pi.est)){
      if(pi.est[i] == 0) counter<-counter+1
    }
    if(counter == 0){
      return(c(NA,NA,NA,NA))
    }
    
    #Compute vector of indices of the prevalences equal to 0 or smaller than pi_min
    indices <- vector() 
    for(i in 1:length(pi.est)){
      if(pi.est[i] < min) indices <- append(indices, i)
    }
    
    c.est <- critpwer(m, pi.est, Sigma, alpha, N)
    
    #Compute estimated c with the minimum number
    c.est_min <- critpwer_min(m, pi.est, Sigma, alpha, N, indices, min)
    
    pwer <- pwerfct(m, pi.true, Sigma, c.est, N)
    pwer_min <- pwerfct(m, pi.true, Sigma, c.est_min, N)
    fwer <- max_fwerfct(m, Sigma, c.est, N)
    fwer_min <- max_fwerfct(m, Sigma, c.est_min, N)
    
    return(c(pwer,pwer_min,fwer,fwer_min))
  }
  cl <- makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  results <- foreach::foreach(seed=1:Nsim) %dopar% {
    set.seed(seed)
    return(f(seed))
  }
  stopCluster(cl)
  return(results)
}


