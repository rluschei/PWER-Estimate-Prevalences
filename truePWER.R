### This R file contains the functions conducting the simulations of the true 
### PWER presented in section 4 of the paper.

library(ggplot2)  
library(MASS)
library(tidyr)
library(doParallel)
no_cores <- detectCores()-2  

### Functions that simulate the values of the true PWER in the various different 
### cases of section 4.

# Prevalences: random and based on independent biomarkers, Estimator: MLE, 
# Variance of the response: unknown
sim1 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    p <- runif(m)
    pi.true <- compute_pi(m,p)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    c.est <- critpwer(m, pi.est, Sigma, alpha, N)
    pwer <- pwerfct(m, pi.true, Sigma, c.est, N)
    return(pwer)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}

# Prevalences: random and based on correlated biomarkers, Estimator: MLE, 
# Variance of the response: unknown
sim2 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    sigma <- randcorr(m)
    rand_vals <- mvrnorm(n=1,mu=rep(0, m),Sigma=sigma)
    p <- pnorm(rand_vals)
    pi.true <- compute_pi(m,p)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    c.est <- critpwer(m, pi.est, Sigma, alpha, N)
    pwer <- pwerfct(m, pi.true, Sigma, c.est, N)
    return(pwer)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}

# Prevalences: equal in all strata, Estimator: MLE, Variance of the response: 
# unknown
sim3 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    pi.true <- matrix(rep(1/(2^m-1),2^m-1), nrow=1, ncol=2^m-1)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    c.est <- critpwer(m, pi.est, Sigma, alpha, N)
    pwer <- pwerfct(m, pi.true, Sigma, c.est, N)
    return(pwer)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}

# Prevalences: 0.5 in one stratum and equal in the others, Estimator: MLE, 
# Variance of the response: unknown
sim4 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    pi.true <- matrix(c(0.5, rep(0.5/(2^m-2),2^m-2)), nrow=1, ncol=2^m-1)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    c.est <- critpwer(m, pi.est, Sigma, alpha, N)
    pwer <- pwerfct(m, pi.true, Sigma, c.est, N)
    return(pwer)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}

# Prevalences: random and based on independent biomarkers, Estimator: Marginal 
# sum estimator, Variance of the response: unknown
sim5 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    p <- runif(m)
    pi.true <- compute_pi(m,p)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    pi.est <- marginal_pi(m, pi.est)
    Sigma <- covmatfct(pi.est, m)
    c.est <- critpwer(m, pi.est, Sigma, alpha, N)
    pwer <- pwerfct(m, pi.true, Sigma, c.est, N)
    return(pwer)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}

# Prevalences: random and based on independent biomarkers, Estimator: MLE, 
# Variance of the response: known
sim6 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    p <- runif(m)
    pi.true <- compute_pi(m,p)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    c.est <- critpwer_kv(m, pi.est, Sigma, alpha)
    pwer <- pwerfct_kv(m, pi.true, Sigma, c.est)
    return(pwer)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}

# Prevalences: random and based on independent biomarkers, Estimator: MLE, 
# Variance of the response: unknown, one single treatment investigated in all
# strata
sim7 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    p <- runif(m)
    pi.true <- compute_pi(m,p)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct_s(pi.est, m)
    c.est <- critpwer(m, pi.est, Sigma, alpha, N)
    pwer <- pwerfct(m, pi.true, Sigma, c.est, N)
    return(pwer)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}


# Execution for m=2-8, N=500
m <- 2:8
df <- sapply(m, sim1, Nsim=10^4, alpha=0.025, N=500)
df <- as.data.frame(df)
names(df) <- paste("m=", as.character(m), sep="")

# Summary of df
summary_df <- data.frame(
  Mean = round(colMeans(df),digits=5),
  SD = round(apply(df, 2, sd),digits=5),
  Min = round(apply(df, 2, min),digits=5),
  Q1 = round(apply(df, 2, quantile, probs = 0.25),digits=5),
  Med = round(apply(df, 2, median),digits=5),
  Q3 = round(apply(df, 2, quantile, probs = 0.75),digits=5),
  Max = round(apply(df, 2, max),digits=5)
)
write.table(summary_df, file = "summary_df.txt", sep="&", quote = F)    

# Plot
df_long <- gather(df, key="m", value="pwer")
plot <- ggplot(df_long, aes(x=m,y=pwer)) + geom_boxplot(outlier.shape=NA) +
  labs(x="",y=expression(PWER(hat(c)))) + coord_cartesian(ylim=c(0.024,0.026)) +
  theme(text = element_text(size = 20))   
ggsave("Plot.pdf", plot, height=8.27, width=11.69)

# Execution for different N, with m=3 fixed
N <- c(25,50,100,150,200,500)
df <- sapply(N, sim1, m=3, Nsim=10^4, alpha=0.025)
df <- as.data.frame(df)
names(df) <- paste("N=", as.character(N), sep="")

# Plot
df_long <- gather(df, key="N", value="pwer")
df_long$N <- factor(df_long$N, levels=c("N=25","N=50","N=100","N=150","N=200",
                                        "N=500"), ordered=T) # Set right order
plot <- ggplot(df_long, aes(x=N,y=pwer)) + geom_boxplot(outlier.shape=NA) +
  labs(x="",y=expression(PWER(hat(c)))) + coord_cartesian(ylim=c(0.02,0.03)) +
  theme(text = element_text(size = 20))   
ggsave("Plot.pdf", plot, height=8.27, width=11.69)




