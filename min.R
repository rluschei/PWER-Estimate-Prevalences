
# Load libraries
source("libraries.R")

# Simulation function
sim <- function(m = 3, Nsim = 10^4, alpha = 0.025, N = 500) {
  cl <- makeCluster(detectCores() - 2)
  registerDoParallel(cl)
  results <- foreach(seed = 1:Nsim, .combine = rbind) %dopar% {
    source("functions.R")
    set.seed(seed)
    p <- runif(m, 0, 0.1)
    tau <- compute_tau(m, p)
    pi.true <- compute_pi(m, tau)
    n <- rmultinom(n = 1, size = N, prob = pi.true)
    
    # Consider only cases where at least one prevalence is 0
    if (!any(n == 0)) return(NA)
    
    pi.est <- n / N
    ALLOC <- allocation(m, n)
    Sigma <- corrmat(m, ALLOC)
    df <- N - sum(ALLOC > 1)
    
    c.est <- critpwer(alpha, m, pi.est, Sigma, df)
    min <- 1 / (2^(m + 1) - 2)
    c.est_min <- critpwer(alpha, m, pi.est, Sigma, df, min)
  
    pwer <- pwer(c.est, m, pi.true, Sigma, df)
    pwer_min <- pwer(c.est_min, m, pi.true, Sigma, df)
    max_swer <- max_swer(c.est, m, pi.true, Sigma, df)
    max_swer_min <- max_swer(c.est_min, m, pi.true, Sigma, df)
    mean_swer <- mean_swer(c.est, m, pi.true, Sigma, df)
    mean_swer_min <- mean_swer(c.est_min, m, pi.true, Sigma, df)
    
    c(pwer, pwer_min, max_swer, max_swer_min, mean_swer, mean_swer_min)
  }
  stopCluster(cl)
  results %>%
    na.omit() %>%
    as.data.frame() %>%
    setNames(c("pwer", "pwer_min", "max_swer", "max_swer_min", "mean_swer", "mean_swer_min"))
}

# Run simulations
df <- sim()

# Create boxplots 
df_long <- df %>%
  pivot_longer(cols = c(pwer, pwer_min, max_swer, max_swer_min, mean_swer, mean_swer_min), 
               names_to = "variable", values_to = "value")

xlabels <- c("unadjusted", "adjusted")
margins <- margin(10, 20, 0, 10)

plot_pwer <- ggplot(df_long %>% filter(variable %in% c("pwer", "pwer_min")), 
                    aes(x = variable, y = value)) + 
  geom_boxplot() + xlab("") +  ylab("true PWER") + scale_x_discrete(labels = xlabels) +
  theme(text = element_text(size = 20), plot.margin = margins) 

plot_max_swer <- ggplot(df_long %>% filter(variable %in% c("max_swer", "max_swer_min")), 
                    aes(x = variable, y = value)) +
  geom_boxplot() + xlab("") + ylab("maximal SWER") + scale_x_discrete(labels = xlabels) +
  theme(text = element_text(size = 20), plot.margin = margins) 

plot_mean_swer <- ggplot(df_long %>% filter(variable %in% c("mean_swer", "mean_swer_min")), 
                         aes(x = variable, y = value)) +
  geom_boxplot() + xlab("") + ylab("mean SWER") + scale_x_discrete(labels = xlabels) +
  theme(text = element_text(size = 20), plot.margin = margins) 

plot <- grid.arrange(plot_pwer, plot_max_swer, plot_mean_swer, ncol = 3)

ggsave(paste0("Plot_min.pdf"), plot, height = 6, width = 15, path = "Results")






