# Load libraries
source("libraries.R")

#' Function for calculation of the PWER
#' 
#' @param c critical value
#' @param m number of populations
#' @param pi prevalence vector
#' @param Sigma correlation matrix
#' @param df degrees of freedom of t-distribution (df = 0 for normal distribution)
#' @param min minimal prevalence
#' 
pwer <- function(c, m, pi, Sigma, df = 0, min = 0){
  pwer <- 0
  below_min <- which(pi < min) 
  if(length(below_min) == 0) fac <- 1 
  else fac <- (1 - length(below_min) * min)/(1 - sum(pi[below_min]))
  for(k in 1:(2^m-1)){
    J <- indexfct(k, m)
    swer_J <- swer(c, Sigma[J,J], df)
    if (k %in% below_min) weight <- min else weight <- pi[k] * fac
    pwer <- pwer + (weight * swer_J)
  }
  pwer
}

#' Function for calculation of the strata-wise FWER
#'
#' @param c critical value
#' @param Sigma correlation matrix for this stratum
#' @param df degrees of freedom of t-distribution
#'
swer <- function(c, Sigma, df = 0){
  if (length(Sigma) == 1) {
    if (df > 0) return(1 - pt(q = c, df = df))
    if (df == 0) return(1 - pnorm(c))
  } else {
    if (df > 0) return(1 - pmvt(upper = rep(c, nrow(Sigma)), sigma = Sigma, df = df)[1])
    if (df == 0) return(1 - pmvnorm(upper = rep(c, nrow(Sigma)), corr = Sigma)[1])
  }
} 
  
#' Function for calculation of the maximal strata-wise FWER
#'
#' @param c critical value
#' @param Sigma correlation matrix for this stratum
#' @param df degrees of freedom of t-distribution
#' 
max_swer <- function(c, m, pi, Sigma, df = 0){
  max <- 0
  for(k in 1:(2^m-1)){
    J <- indexfct(k, m)
    swer <- swer(c, Sigma[J,J], df)
    if(swer > max & pi[k] > 0) max <- swer
  }
  max
}

#' Function for calculation of the mean strata-wise FWER
#'
#' @param c critical value
#' @param Sigma correlation matrix for this stratum
#' @param df degrees of freedom of t-distribution
#' 
mean_swer <- function(c, m, pi, Sigma, df = 0){
  mean <- 0
  for(k in 1:(2^m-1)){
    J <- indexfct(k, m)
    if(pi[k] > 0) mean <- mean + swer(c, Sigma[J,J], df)
  }
  mean/(2^m-1)
}

#' Returns the index set J belonging to the k-th stratum (containing the relevant 
#' treatments for this stratum)
#' 
#' @param k number(s) of stratum (between 1 and 2^m-1). Can be a scalar or a vector.
#' @param m number of populations
#' 
indexfct <- function(k, m){
  # Define a helper function for individual k
  index_single <- function(k_val, m){
    v <- binvector(k_val, m)
    set <- vector()
    for(i in 1:(length(v))) {
      if(v[i] == 1) set <- append(set, i)
    }  
    set
  }
  
  # Apply the function element-wise if k is a vector
  if (length(k) > 1) {
    lapply(k, function(k_val) index_single(k_val, m))
  } else {
    index_single(k, m)
  }
}


#' Conversion of the natural number k into a binary vector with d digits
#' 
#' @param k a natural number
#' @param d the number of digits
#' 
binvector <- function(k, d) {  
  vec = rev(as.numeric(intToBits(k)))
  vec[-(1:(length(vec) - d))]
}

#' Computes the correlation matrix
#' 
#' @param m number of biomarkers
#' @param ALLOC allocation matrix
#'
corrmat <- function(m, ALLOC, VAR = NULL){
  Sigma <- matrix(rep(1, m*m), nr = m)
  for(i in 1:m){
    for(j in 1:m){
      if(i != j) Sigma[i,j] <- corr(i, j, m, ALLOC, VAR)
    }
  } 
  Sigma
}

#' Calculation of the correlation between two given test statistics
#' 
#' @param i index of first test statistic
#' @param j index of second test statistic
#' @param m total number of biomarkers
#' @param ALLOC allocation matrix with 2^m-1 rows and m+1 columns. In entry ALLOC[i,j] 
#'        it contains the sample size for treatment j in stratum i. The last 
#'        column is for the control.
#'
corr <- function(i, j, m, ALLOC, VAR = NULL) {
  single <- ncol(ALLOC) == 2 # true in single treatment case
  
  n_iC <- n_i(i, m, ALLOC, "C", single)
  n_jC <- n_i(j, m, ALLOC, "C", single) 
  n_iT <- n_i(i, m, ALLOC, "T", single) 
  n_jT <- n_i(j, m, ALLOC, "T", single)
  
  if (n_iC == 0 | n_jC == 0 | n_iT == 0 | n_jT == 0) return(0)
  
  if(is.null(VAR)) V_i <- 1/n_iT + 1/n_iC else V_i <- V(i, m, ALLOC, VAR, single)
  if(is.null(VAR)) V_j <- 1/n_jT + 1/n_jC else V_j <- V(j, m, ALLOC, VAR, single)
  
  res <- 0
  for(k in 1:(2^m - 1)) {
    if(i %in% indexfct(k, m) & j %in% indexfct(k, m)) {
      if (single) {
          if(is.null(VAR)) res <- res + ALLOC[k, 1] / (n_iT * n_jT) + ALLOC[k, 2] / (n_iC * n_jC)
          else res <- res + ALLOC[k, 1] * VAR[k, 1] / (n_iT * n_jT) + ALLOC[k, 2] * VAR[k, 2] / (n_iC * n_jC)
      } else {
          if(is.null(VAR)) res <- res + ALLOC[k, m+1] / (n_iC * n_jC)
          else res <- res + ALLOC[k, m+1] * VAR[k, m+1] / (n_iC * n_jC)
      }
    }
  }
  res / sqrt(V_i * V_j)
}

# Returns the sample size in population i for treatment tmt
n_i <- function(i, m, ALLOC, tmt, single){
  sum <- 0
  for(k in 1:(2^m - 1)) {
    if(i %in% indexfct(k, m)) {
      if(tmt == "C") sum <- sum + (if (single) ALLOC[k, 2] else ALLOC[k, m+1])
      else sum <- sum + (if (single) ALLOC[k, 1] else ALLOC[k, i])
    }
  }
  sum
}


#' Calculation of the critical value for PWER-control
#' 
#' @param m number of populations
#' @param pi prevalence vector
#' @param ALLOC allocation matrix
#' @param alpha significance level
#' @param distr distribution of test statistics ("z" or "t")
#' @param min minimum prevalence 
#' 
critpwer <- function(alpha, m, pi, Sigma, df = 0, min = 0){
  f <- function(c){pwer(c, m, pi, Sigma, df, min) - alpha}
  uniroot(f, interval = c(0,20))$root
}

#' Calculation of the critical value for FWER-control
#' 
#' @param m number of populations
#' @param ALLOC allocation matrix
#' @param alpha significance level
#' @param distr distribution of test statistics ("z" or "t")
#' 
critfwer <- function(alpha, m, Sigma, df = 0){
  f <- function(c){swer(c, Sigma, df) - alpha}
  uniroot(f, interval = c(0,20))$root
}


#' Calculates the prevalences of the relevant strata
#' 
#' @param m number of biomarkers
#' @param tau prevalence vector, including the stratum where no biomarkers are 
#' present in first element
#' 
compute_pi <- function(m, tau){
  tau[2:(2^m)]/(1-tau[1])
}

#' Calculates the prevalences from the biomarker expression probabilities 
#' (including the stratum where none of the biomarkers are present)
#' 
#' @param m number of biomarkers
#' @param p vector of biomarker expression probabilities
#' 
compute_tau <- function(m, p){
  tau <- rep(1, 2^m)
  for(k in 1:(2^m)){
    bin <- binvector(k-1, m)
    for(i in 1:m){
      if(bin[i] == 1){
        tau[k] <- tau[k] * p[i]
      }
      else{
        tau[k] <- tau[k] * (1 - p[i])
      }
    }
  }
  tau
}


#' Computes the vector of marginal sum estimators
#' 
#' @param m number of biomarkers
#' @param tau prevalence vector (including the stratum where none of the 
#' biomarkers are present)  
#' 
compute_marg_est <- function(m, tau){
  p.est <- rep(0, m)
  for(i in 1:m){
    for(k in 1:(2^m-1)){
      J <- indexfct(k, m)
      if(i %in% J) p.est[i] <- p.est[i] + tau[k+1]
    }
  }
  tau.marg <- compute_tau(m, p.est)
  compute_pi(m, tau.marg)
}

#' Creates a matrix containing allocation numbers for the given allocation
#' strategy. The rows of the matrix correspond to the strata (according to the 
#' indexfct) and the columns to the treatments. The last column is for the control.
#' 
#' @param m number of biomarkers
#' @param n vector of strata-wise sample sizes
#' @param strategy allocation strategy ("equal" for equal allocation probabilities within the strata, 
#'        "random" for random allocation numbers based on equal allocation probabilities,
#'        "random_probs" for random allocation numbers based on random allocation probabilities)
#' @param single indicating if treatments are equal (TRUE) or pairwise different (FALSE)
#'        
allocation <- function(m, n, strategy = "equal", single = F) { 
  nrow <- 2^m - 1
  ncol <- if (single) 2 else m + 1
  ALLOC <- matrix(0, nrow = nrow, ncol = ncol)
  
  for (i in 1:nrow) {
    tmts <- if (single) c(1, 2) else c(indexfct(i, m), m + 1)
    probs <- rep(1/length(tmts), length(tmts))
    
    if (strategy == "random_probs") {
      probs <- runif(length(tmts))
      probs <- probs / sum(probs)
    }
    
    if (strategy == "equal") {
      alloc <- floor(n[i] * probs)
      remainder <- n[i] - sum(alloc)
      ALLOC[i, tmts] <- alloc + rmultinom(1, remainder, probs)
    } else {
      ALLOC[i, tmts] <- rmultinom(1, n[i], probs)
    }
  }
  ALLOC
}

#' Simulation function that computes the values of the true PWER, the maximal SWER 
#' or the mean SWER in the various different cases from section 4.
#'
#' @param m number of biomarkers
#' @param N sample size
#' @param Nsim number of simulation runs
#' @param alpha significance level for PWER-control
#' @param indep_bio indicates if biomarkers are independent (TRUE) or not (FALSE) 
#' @param prev_mode indicates if the prevalences should be randomly generated ("random"), 
#' fixed and equal ("fixed_equal") or fixed with one big value ("fixed_0.5")
#' @param distr distribution of test statistics ("z" or "t")
#' @param strategy allocation strategy (see allocation function)
#' @param single indicating if treatments are equal (TRUE) or pairwise different (FALSE)
#' @param est estimator ("MLE" or "marg_est")
#' @param error_rate error rate ("pwer", "max_swer" or "mean_swer")
#'
sim <- function(m, N, Nsim = 10^4, alpha, indep_bio, prev_mode, distr, strategy,
                single, est, error_rate = "pwer", het_var = F){
  cl <- makeCluster(detectCores() - 2)
  registerDoParallel(cl)
  results <- foreach(seed = 1:Nsim, .combine = c) %dopar% {
    source("functions.R")
    set.seed(seed)
    if(prev_mode == "fixed_equal") tau <- matrix(c(0, rep(1/(2^m - 1), 2^m - 1)), 
                                                 nrow = 1, ncol = 2^m)
    if(prev_mode == "fixed_0.5") tau <- matrix(c(0, 0.5, rep(0.5/(2^m - 2), 2^m - 2)), 
                                               nrow = 1, ncol = 2^m)
    if(prev_mode == "random"){
      if(indep_bio) p <- runif(m)
      else {
        sigma <- randcorr(m)
        rand_vals <- mvrnorm(n = 1, mu = rep(0,m), Sigma = sigma)
        p <- pnorm(rand_vals)
      }
      tau <- compute_tau(m, p)
    }
    pi.true <- compute_pi(m, tau)
    n <- rmultinom(n = 1, size = N, prob = tau)
    if(n[1] == N) return(NA) # No screened patients
    if(est == "MLE") pi.est <- n[-1]/(N - n[1]) 
    if(est == "marg_est") {tau.est <- n/N; pi.est <- compute_marg_est(m, tau.est)}
    ALLOC <- allocation(m, n[-1], strategy, single)
    if(distr == "t") df <- N - sum(ALLOC > 1) else df <- 0
    if(het_var) VAR <- gen_var(m, single) else VAR <- NULL
    Sigma <- corrmat(m, ALLOC, VAR)
    c.est <- critpwer(alpha, m, pi.est, Sigma, df)
    if(error_rate == "pwer") return(pwer(c.est, m, pi.true, Sigma, df))
    if(error_rate == "max_swer") return(max_swer(c.est, m, pi.true, Sigma, df))
    if(error_rate == "mean_swer") return(mean_swer(c.est, m, pi.true, Sigma, df))
  }
  stopCluster(cl)
  results
}

# Summary table function
summary_tbl <- function(df, error_rate, case_no){
  stats <- c("Mean" = "{mean}", "SD" = "{sd}", "Min" = "{min}", "Q1" = "{p25}",
             "Median" = "{median}", "Q3" = "{p75}", "Max" = "{max}")
  tbl <- 
    imap(
      stats,
      ~df %>%
        tbl_summary(type = list(everything() ~ "continuous"), 
                    statistic = everything()~.x, digits = all_continuous() ~ 5, 
                    missing = "no")  %>%
        modify_header(all_stat_cols() ~ stringr::str_glue("**{.y}**"))
    ) %>%
    tbl_merge(tab_spanner = FALSE) %>% modify_footnote(~NA) %>% 
    modify_header(label ~ ifelse(case_no == 11, "**N**", "**m**")) %>% as_gt() %>% 
    gtsave(paste0("Summary_", error_rate, "_case", case_no, ".pdf"), path = "Results")
}

#' Generates a sample of observations under the null hypothesis 
#'
#' @param m number of biomarkers
#' @param ALLOC allocation matrix
#' @param VAR matrix of strata-treatment wise variances
#' 
gen_sample <- function(m, ALLOC, VAR) {
  # Convert ALLOC and VAR matrices to data frames for easier manipulation and 
  # add Stratum as a row number column 
  ALLOC_df <- as.data.frame(ALLOC) %>%
    set_names(c(paste0("T", 1:(ncol(.) - 1)), "C")) %>%
    mutate(Stratum = row_number())
  
  VAR_df <- as.data.frame(VAR) %>% 
    set_names(c(paste0("T", 1:(ncol(.) - 1)), "C")) %>%
    mutate(Stratum = row_number())
  
  # Convert both data frames to long format
  ALLOC_long <- ALLOC_df %>% pivot_longer(-Stratum, names_to = "Treatment", values_to = "Allocated")
  VAR_long <- VAR_df %>% pivot_longer(-Stratum, names_to = "Treatment", values_to = "Variance")
  
  # Join the two data frames by Stratum and Treatment
  combined_df <- left_join(ALLOC_long, VAR_long, by = c("Stratum", "Treatment"))
  
  # Determine the population indicators for each stratum
  population_indicators <- data.frame(Stratum = unique(combined_df$Stratum))
  for (i in 1:m) {
    population_indicators[[paste0("P", i)]] <- 0
  }
  
  for (k in unique(population_indicators$Stratum)) {
    indices <- indexfct(k, m)
    population_indicators[population_indicators$Stratum == k, paste0("P", indices)] <- 1
  }
  
  # Combine the indicators with the original data
  combined_df <- combined_df %>%
    left_join(population_indicators, by = "Stratum")
  
  # Generate the random samples, with number of samples based on ALLOC (Allocated)
  sample <- combined_df %>%
    rowwise() %>% # Work on each row
    mutate(Response = list(rnorm(Allocated, mean = 0, sd = sqrt(Variance)))) %>% # Generate multiple random samples
    unnest(Response) %>% # Unnest the list of responses into individual rows
    select(Stratum, Treatment, Response, everything(), -Allocated, -Variance)
}

#' Generates random strata-treatment-specific residual variances 
#' 
#' @param m number of populations
#' @param single indicates if treatments are equal (TRUE) or pairwise different (FALSE)
#' 
gen_var <- function(m, single){
  nrow = 2^m-1
  if (single) ncol <- 2 else ncol <- m + 1
  VAR <- matrix(NA, nrow = nrow, ncol = ncol)
  for (k in 1:nrow) {
    if(single) columns <- c(1,2)
    else columns <- c(indexfct(k, m), ncol)
    VAR[k, columns] <- runif(length(columns), min = 0, max = 1)
  }
  VAR
}

# Compute estimated strata-treatment-wise variances from a sample
var.est <- function(sample) {
  sample %>%
    group_by(Stratum, Treatment) %>%
    summarise(Variance = var(Response, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Treatment, values_from = Variance, names_sort = TRUE) %>%
    select(Stratum, everything(), -starts_with("C"), "C") %>% # Move column "C" to the end
    column_to_rownames(var = "Stratum") %>%
    as.matrix()
}

#' Compute degrees of freedom according to Satterthwaite
#' 
#' @param i population index
#' @param m number of populations
#' @param sample sample
#' @param single indicates if treatments are equal (TRUE) or pairwise different (FALSE)
#' 
df <- function(i, m, sample, single){
  sample_i <- sample[sample[[paste0("P", i)]] == 1,]
  sample_iC <- sample_i[sample_i$Treatment == "C", ]$Response
  if(single) sample_iT <- sample_i[sample_i$Treatment == "T1", ]$Response
  else sample_iT <- sample_i[sample_i$Treatment == paste0("T", i), ]$Response
  
  n_iT <- length(sample_iT)
  n_iC <- length(sample_iC)
  
  num <- (var(sample_iT)/n_iT + var(sample_iC)/length(sample_iC))^2
  denom <- (var(sample_iT)/n_iT)^2/(n_iT - 1) + 
    (var(sample_iC)/n_iC)^2/(n_iC - 1)
  
  floor(num/denom)
}

#' Calculation of the critical values for PWER-control for the unknown, 
#' heterogeneous variances case
#' 
#' @param alpha significance level
#' @param m number of populations
#' @param pi prevalence vector
#' @param Sigma.est estimated correlation matrix
#' @param sample sample
#' @param single indicates if treatments are equal (TRUE) or pairwise different (FALSE)
#' 
critpwer_sw <- function(alpha, m, pi, Sigma.est, sample, single){
  crit <- rep(0, m)
  for(i in 1:m){
    df_i = df(i, m, sample, single)
    f <- function(c){pwer(c, m, pi, Sigma.est, df_i) - alpha}
    crit[i] <- uniroot(f, interval = c(0,20))$root
  }
  crit
}


#' Computes the variance of the test statistics
#' 
#' @param i index of test statistic
#' @param m number of populations
#' @param ALLOC allocation matrix
#' @param VAR matrix of strata-treatment-wise variances
#' 
V <- function(i, m, ALLOC, VAR, single){ 
  n_iC <- n_i(i, m, ALLOC, "C", single)
  n_iT <- n_i(i, m, ALLOC, "T", single) 
  res <- 0
  for(k in 1:(2^m - 1)) {
    if(i %in% indexfct(k, m)) {
      if(single) res <- res + ALLOC[k, 1] * VAR[k, 1] / n_iT^2 + ALLOC[k, 2] * VAR[k, 2] / n_iC^2
      else res <- res + ALLOC[k, i] * VAR[k, i] / n_iT^2 + ALLOC[k, m+1] * VAR[k, m+1] / n_iC^2
    }
  }
  res
}

#' Computes the test statistics from a sample
#' 
#' @param m number of populations
#' @param sample sample
#' @param ALLOC allocation matrix
#' @param VAR matrix of strata-treatment-wise variances
#' 
test <- function(m, sample, ALLOC, VAR){
  single = single <- ncol(ALLOC) == 2 # true in single treatment case
  retv <- rep(0, m)
  for (i in 1:m){
    sample_i <- sample[sample[[paste0("P", i)]] == 1,]
    sample_iC <- sample_i[sample_i$Treatment == "C", ]$Response
    if(single) sample_iT <- sample_i[sample_i$Treatment == "T1", ]$Response
    else sample_iT <- sample_i[sample_i$Treatment == paste0("T", i), ]$Response
      
    retv[i] <- (mean(sample_iT) - mean(sample_iC))/sqrt(V(i, m, ALLOC, VAR, single))
  }
  retv
}

# Computes the PWER approximation from the test results
pwer_sim <- function(m, pi, rejections, Nsim){
  pwer <- 0
  for(k in 1:(2^m-1)){
    J <- indexfct(k,m)
    if(length(J) > 1) pwer <- pwer + pi[k] * sum(rowSums(rejections[, indexfct(k,m)] > 0) > 0)
    else pwer <- pwer + pi[k] * sum(rejections[, indexfct(k,m)] > 0)
  }
  pwer/Nsim
}


