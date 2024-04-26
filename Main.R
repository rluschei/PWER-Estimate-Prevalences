### This R file contains the functions conducting the simulations presented
### in section 4.

library(doParallel)
no_cores <- detectCores()-2  
library(foreach)
library(purrr)

# Simulation function that computes the values of the PWER and the maximal FWER 
#in the various different cases from section 4.
sim <- function(m, N, Nsim=10^4, alpha=0.025, value="pwer", indep_bio=T, 
                prevs="rand", distr="t", alloc="equal", case="d", est="MLE"){
  cl <- makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  results <- foreach::foreach(seed=1:Nsim) %dopar% {
    f <- function(seed, m, alpha, N, value, indep_bio, prevs, distr, alloc, 
                  case, est){
      source("Functions.R")
      if(prevs=="fixed_eq") tau <- matrix(c(0,rep(1/(2^m-1),2^m-1)), nrow=1, ncol=2^m)
      if(prevs=="fixed_0.5") tau <- matrix(c(0,0.5, rep(0.5/(2^m-2),2^m-2)), nrow=1, ncol=2^m)
      if(prevs=="rand"){
        if(indep_bio) p <- runif(m)
        else{
           sigma <- randcorr(m)
           rand_vals <- mvrnorm(n=1,mu=rep(0,m),Sigma=sigma)
          p <- pnorm(rand_vals)
        }
        tau <- compute_tau(m,p)
      }
      pi.true <- compute_pi(m,tau)
      n <- rmultinom(n=1, size=N, prob=tau)
      if(n[1]==N) return(c(NA,NA)) # No screened patients
      if(distr=="t"&(N-n[1]<=2^m-1)) return(c(NA,NA)) #Not enough patients for positive df
      if(est=="MLE") pi.est <- n[-1]/(N-n[1]) 
      if(est=="marg") {tau.est <- n/N; pi.est <- marginal_pi(m,tau.est)}
      ALLOC <- allocation(m,n[-1],alloc,case)
      c.est <- critpwer(m, pi.est, ALLOC, alpha, distr)
      pwer <- pwer(m, pi.true, ALLOC, c.est, distr)
      fwer <- max_fwer(m, ALLOC, c.est, distr)
      return(c(pwer,fwer))
    }
    set.seed(seed)
    return(f(seed,m,alpha,N,value,indep_bio,prevs,distr,alloc,case,est))
  }
  stopCluster(cl)
  return(unlist(results))
}


### Running the simulations

source("Exportfunctions.R")

### True PWER and maximal strata-wise FWER

# Calls the sim function with the given arguments and exports the results
wrap <- function(m=2:8, N=500, value="pwer", indep_bio=T, prevs="rand", distr="t", alloc="equal", 
              case="d", est="MLE", i){
  if(length(m)>1) df_temp <- sapply(m, sim, N=N, value=value, indep_bio=indep_bio, 
                                    prevs=prevs, distr=distr, alloc=alloc, 
                                    case=case, est=est)
  if(length(N)>1) df_temp <- sapply(N, sim, m=m, value=value, indep_bio=indep_bio, 
                                     prevs=prevs, distr=distr, alloc=alloc, 
                                     case=case, est=est)
  export(m,N,i,df_temp)
  return(df_temp)
} 

# List of the different cases
args <- list(
  list(i=1),
  list(i=2, indep_bio = FALSE),
  list(i=3, prevs = "fixed_eq"),
  list(i=4, prevs = "fixed_0.5"),
  list(i=5, est = "marg"),
  list(i=6, distr = "z"),
  list(i=7, case = "s"),
  list(i=8, alloc = "randomp"),
  list(i=9, m=3, N=c(25,50,100,150,200,500)),
  list(i=10, m=4, N=c(25,50,100,150,200,500))
)

df <- map(args, ~do.call(wrap, .x))



#### Minimum number 
sim_min <- function(m, Nsim=10^4, alpha=0.025, N=500){
  cl <- makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  results <- foreach::foreach(seed=1:Nsim) %dopar% {
    f <- function(seed){
      source("Functions.R")
      p <- runif(m,0,0.1)
      tau <- compute_tau(m,p)
      pi.true <- compute_pi(m,tau)
      n <- rmultinom(n=1, size=N, prob=pi.true)
      if(n[1]==N) return(c(NA,NA,NA,NA))
      if(N-n[1]<=2^m-1) return(c(NA,NA,NA,NA)) 
      pi.est <- n/N
      
      # test if one prevalence is 0
      counter = 0
      for(i in 1:length(pi.est)){
        if(pi.est[i] == 0) counter<-counter+1
      }
      if(counter == 0){
        return(c(NA,NA,NA,NA))
      }
      
      ALLOC <- allocation(m,n)
      c.est <- critpwer(m, pi.est, ALLOC, alpha)
      
      #Compute estimated c with the minimum number
      min <- 1/(2^(m+1)-2)
      c.est_min <- critpwer(m, pi.est, ALLOC, alpha, min=min)
      
      pwer <- pwer(m, pi.true, ALLOC, c.est)
      pwer_min <- pwer(m, pi.true, ALLOC, c.est_min)
      fwer <- max_fwer(m, ALLOC,c.est)
      fwer_min <- max_fwer(m, ALLOC, c.est_min)
      
      return(c(pwer,pwer_min,fwer,fwer_min))
    }
    set.seed(seed)
    return(f(seed))
  }
  stopCluster(cl)
  results <- as.data.frame(results)
  boxplots_min(results, m)
  summary_table_min(results, m)
  return(results)
}


df_min <- sapply(2:8, sim_min)
  

 
  











