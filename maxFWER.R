### This R file contains the functions conducting the simulations of the  
### maximal strata-wise FWER presented in section 4 of the paper.

library(ggplot2)  
library(MASS)
library(tidyr)
library(doParallel)
no_cores <- detectCores()-2  

# Simulates the strata-wise FWER with the estimated rejection boundary
sim1 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    p <- runif(m)
    pi.true <- compute_pi(m,p) 
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    c <- critpwer(m, pi.est, Sigma, alpha, N)
    max <- max_fwerfct(m, Sigma, c, N)
    return(max)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}

# Simulates the strata-wise FWER with the true rejection boundary
sim2 <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    p <- runif(m)
    pi.true <- compute_pi(m,p) 
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    c <- critpwer(m, pi.true, Sigma, alpha, N)
    max <- max_fwerfct(m, Sigma, c, N)
    return(max)
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}


m <- 2:4
df1 <- sapply(m, sim1, Nsim=10^2, alpha=0.025, N=500)
df1 <- as.data.frame(df1)
names(df1) <- paste("m=", as.character(m), sep="")

df2 <- sapply(m, sim2, Nsim=10^2, alpha=0.025, N=500)
df2 <- as.data.frame(df2)
names(df2) <- paste("m=", as.character(m), sep="")

diff <- df1-df2

# Plot
df_long <- gather(df1, key="m", value="fwer")
plot <- ggplot(df_long, aes(x=m,y=fwer)) + geom_boxplot(outlier.shape=NA) +
  labs(x="",y=expression(max[J %subseteq% I]~FWER[J](hat(c)))) + 
  coord_cartesian(ylim=c(0.024,0.08)) +
  theme(text = element_text(size = 20))   
ggsave("Plot.pdf", plot, height=8.27, width=11.69)

# Plot
df_long <- gather(df2, key="m", value="fwer")
plot <- ggplot(df_long, aes(x=m,y=fwer)) + geom_boxplot(outlier.shape=NA) +
  labs(x="",y=expression(abs(max[J %subseteq% I]~(FWER[J](hat(c))-FWER[J](c))))) + 
  coord_cartesian(ylim=c(-0.05,0.05)) +
  theme(text = element_text(size = 20))   
ggsave("Plot.pdf", plot, height=8.27, width=11.69)

