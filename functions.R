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
#' @param k number of stratum (between 1 and 2^m-1)
#' @param m number of populations
#' 
indexfct <- function(k, m){
  v <- binvector(k, m)
  set <- vector()
  for(i in 1:(length(v))) {
    if(v[i] == 1) set <- append(set, i)
  }  
  set
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
corrmat <- function(m, ALLOC){
  Sigma <- matrix(rep(1, m*m), nr = m)
  for(i in 1:m){
    for(j in 1:m){
      if(i != j) Sigma[i,j] <- corr(i, j, m, ALLOC)
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
corr <- function(i, j, m, ALLOC) {
  n_iC <- 0; n_jC <- 0; n_iT <- 0; n_jT <- 0
  single <- ncol(ALLOC) == 2 # true in single treatment case
  for(k in 1:(2^m - 1)) {
    if(i %in% indexfct(k, m)) {
      n_iC <- n_iC + (if (single) ALLOC[k, 2] else ALLOC[k, m+1])
      n_iT <- n_iT + (if (single) ALLOC[k, 1] else ALLOC[k, i])
    }
    if(j %in% indexfct(k, m)) {
      n_jC <- n_jC + (if (single) ALLOC[k, 2] else ALLOC[k, m+1])
      n_jT <- n_jT + (if (single) ALLOC[k, 1] else ALLOC[k, j])
    }
  }
  if (n_iC == 0 | n_jC == 0 | n_iT == 0 | n_jT == 0) return(0)

  H_i <- 1/n_iT + 1/n_iC
  H_j <- 1/n_jT + 1/n_jC

  res <- 0
  for(k in 1:(2^m - 1)) {
    if(i %in% indexfct(k, m) & j %in% indexfct(k, m)) {
      if (single) {
        res <- res + ALLOC[k, 1] / (n_iT * n_jT) + ALLOC[k, 2] / (n_iC * n_jC)
      } else {
        res <- res + ALLOC[k, m+1] / (n_iC * n_jC)
      }
    }
  }
  res / sqrt(H_i * H_j)
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
#' @param tmt_comp indicating if treatments are equal ("equal) or pairwise different ("different")
#' @param strategy allocation strategy ("equal" for equal allocation probabilities within the strata, 
#'        "random" for random allocation numbers based on equal allocation probabilities,
#'        "random_probs" for random allocation numbers based on random allocation probabilities)
#'        
allocation <- function(m, n, strategy = "equal", tmt_comp = "different") {
  nrow <- 2^m - 1
  ncol <- if (tmt_comp == "different") m + 1 else 2
  ALLOC <- matrix(0, nrow = nrow, ncol = ncol)
  
  for (i in 1:nrow) {
    tmts <- if (tmt_comp == "different") c(indexfct(i, m), m + 1) else c(1, 2)
    probs <- rep(1/length(tmts), length(tmts))
    
    if (strategy == "random_probs") {
      probs <- runif(length(tmts))
      probs <- probs / sum(probs)
    }
    
    if (strategy == "equal") {
      ALLOC[i, tmts] <- n[i] * probs
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
#' @param tmt_comp indicating if treatments are equal ("equal) or pairwise different ("different")
#' @param est estimator ("MLE" or "marg_est")
#' @param error_rate error rate ("pwer", "max_swer" or "mean_swer")
#'
sim <- function(m, N, Nsim = 10^4, alpha, indep_bio, prev_mode, distr, strategy,
                tmt_comp, est, error_rate = "pwer"){
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
    ALLOC <- allocation(m, n[-1], strategy, tmt_comp)
    if(distr == "t") df <- N - sum(ALLOC > 1) else df <- 0
    Sigma <- corrmat(m, ALLOC)
    
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
    modify_header(label ~ ifelse(case_no == 10, "**N**", "**m**")) %>% as_gt() %>% 
    gtsave(paste0("Summary_", error_rate, "_case", case_no, ".pdf"), path = "Results")
}


