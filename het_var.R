
# Load libraries
source("libraries.R")

# Simulation function
sim <- function(m = 2, Nsim = 10^2, alpha = 0.025, N = 500) {
  cl <- makeCluster(detectCores() - 2)
  registerDoParallel(cl)
  results <- foreach(seed = 1:Nsim, .combine = rbind) %dopar% {
    source("functions.R")
    set.seed(seed)
    
    # Generate true prevalences
    p <- runif(m)
    tau <- compute_tau(m, p)
    pi.true <- compute_pi(m, tau)
    
    # Generate strata-wise sample sizes and compute MLE of prevalences
    n <- rmultinom(n = 1, size = N, prob = pi.true)
    pi.est <- n / N

    # Generate treatment allocation matrix
    single <- F
    ALLOC <- allocation(m, n, single = single)

    # Generate strata-treatment specific variances
    VAR <- gen_var(m, single)

    # Define matrix to store the test decisions
    rejections <- matrix(0, nrow = Nsim, ncol = m)

    # Approximate the true PWER
    for(i in 1:Nsim){
      #Generate sample and variance estimations
      sample <- gen_sample(m, ALLOC, VAR)
      VAR.est <- var.est(sample)

      # Return NA if not all strata-treatment-wise variances can be estimated
      if(sum(!is.na(VAR)) !=  sum(!is.na(VAR.est))) return(NA)

      # Compute correlation matrix
      Sigma.est <- corrmat(m, ALLOC, VAR.est)

      # Compute critical values from Satterthwaite approximation
      c.est <- critpwer_sw(alpha, m, pi.est, Sigma.est, sample, single)

      # Compute test statistics
      t_star <- test(m, sample, ALLOC, VAR.est)

      # Store test decisions
      rejections[i,] <- t_star > c.est
    }
    pwer <- pwer_sim(m, pi.true, rejections, Nsim)
  }
  stopCluster(cl)
  results
}

# Run simulations
df <- as.data.frame(sim())

# Create boxplot
plot <- ggplot(df, aes(y = V1)) + 
  geom_boxplot() +
  labs(x = "", y = expression(PWER(hat(c)))) + 
  theme(text = element_text(size = 20),
        axis.title.x = element_text(margin = margin(t = 10)), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(paste0("Plot_pwer_approx.pdf"), plot, height = 6, width = 8, path = "Results")






