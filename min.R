library(ggplot2)  
library(MASS)
library(tidyr)
library(doParallel)
no_cores <- detectCores()-2
library(gridExtra)

sim <- function(m, Nsim, alpha, N){
  f <- function(k){
    source("Functions.R")
    p <- runif(m,0,0.2)
    pi.true <- compute_pi(m,p)
    pi.est <- rmultinom(n=1, size=N, prob=pi.true)/N
    Sigma <- covmatfct(pi.est, m)
    
    indices <- vector() #Vector for the indices of the sample sizes equal to 0
    for(i in 1:length(pi.est)){
      if(pi.est[i] == 0) indices <- append(indices, i)
    }
    if(length(indices) == 0){
      return(c(NA,NA,NA,NA))
    }
    
    c.est <- critpwer(m, pi.est, Sigma, alpha, N)
    
    #Compute estimated c with the minimum number
    c.est_min <- critpwer_min(m, pi.est, Sigma, alpha, N, indices)
    
    pwer <- pwerfct(m, pi.true, Sigma, c.est, N)
    pwer_min <- pwerfct(m, pi.true, Sigma, c.est_min, N)
    fwer <- max_fwerfct(m, Sigma, c.est, N)
    fwer_min <- max_fwerfct(m, Sigma, c.est_min, N)
    
    return(c(pwer,pwer_min,fwer,fwer_min))
  }
  cl <- makeCluster(no_cores)
  results <- parSapply(cl, 1:Nsim, f)
  stopCluster(cl)
  return(results)
}

df <- sim(m=3, Nsim=10^4, alpha=0.025, N=500)
df <- as.data.frame(df)
df$group <- c("PWER","PWER_min","FWER","FWER_min")
df_sub1 <- df[1:2,]
df_long1 <- gather(df_sub1, key = "variable", value = "value", -group)
df_sub2 <- df[3:4,]
df_long2 <- gather(df_sub2, key = "variable", value = "value", -group)

plot1 <- ggplot(df_long1, aes(x = group, y = value)) + 
  geom_boxplot(outlier.shape=NA) + xlab("") + ylab("") + 
  scale_x_discrete(labels=c(expression(PWER(hat(c))), expression(PWER(hat(c)[min])))) +
  theme(text = element_text(size = 20)) 
plot2 <- ggplot(df_long2, aes(x = group, y = value)) +
  geom_boxplot(outlier.shape=NA) + xlab("") + ylab("") +
scale_x_discrete(labels=c(expression(max[J %subseteq% I]~FWER[J](hat(c))), expression(max[J %subseteq% I]~FWER[J](hat(c)[min])))) +
  theme(text = element_text(size = 20)) 
plot <- grid.arrange(plot1, plot2, ncol = 2)
ggsave("Plot.pdf", plot, height=6, width=11.69)


