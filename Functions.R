### This R file contains the functions that are called in the files truePWER.R,
### maxSWER.R and min.R. 


library(multcomp)
library(randcorr)

# Function for calculation of the PWER
pwerfct <- function(m, piv, Sigma, c, N){
  pwer = 0
  for(k in 1:(2^m-1)){
    J <- indexfct(k, m)
    card_J <- length(J)
    Sigma_J <- subcovmatfct(J, card_J, Sigma)
    n <- 2^m-1
    if(card_J==1){
      pwer <- pwer + (piv[k] * (1-pt(q=c, df=N-n)))
    }
    else{
      pwer <- pwer + (piv[k] * (1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1])) 
    }
  } 
  return (pwer)
}

# Function for calculation of the PWER in case of known variances of the 
#response
pwerfct_kv <- function(m, piv, Sigma, c){
  pwer = 0
  for(k in 1:(2^m-1)){
    J <- indexfct(k, m)
    card_J <- length(J)
    Sigma_J <- subcovmatfct(J, card_J, Sigma)
    if(card_J==1){
      pwer <- pwer + (piv[k] * (1-pnorm(c)))
    }
    else{
      pwer <- pwer + (piv[k] * (1-pmvnorm(upper=rep(c,card_J), corr=Sigma_J)[1])) 
    }
  } 
  return (pwer)
}

# Function for calculation of the maximal strata-wise FWER
max_fwerfct <- function(m, Sigma, c, N){
  max = 0
  for(k in 1:(2^m-1)){
    J <- indexfct(k, m)
    card_J <- length(J)
    Sigma_J <- subcovmatfct(J, card_J, Sigma)
    n <- 2^m-1
    if(card_J==1){
      if(1-pt(q=c, df=N-n)>max) max <- 1-pt(q=c, df=N-n)
    }
    else{
      if(1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1]>max){
        max <- 1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1] 
      }
    }
  } 
  return (max)
}

# Function for calculation of the maximal strata-wise FWER
max_fwerfct_kv <- function(m, Sigma, c){
  max = 0
  for(k in 1:(2^m-1)){
    J <- indexfct(k, m)
    card_J <- length(J)
    Sigma_J <- subcovmatfct(J, card_J, Sigma)
    n <- 2^m-1
    if(card_J==1){
      if(1-pnorm(c)>max) max <- 1-pnorm(c)
    }
    else{
      if(1-pmvnorm(upper=rep(c,card_J), corr=Sigma_J)[1]>max){
        max <- 1-pmvnorm(upper=rep(c,card_J), corr=Sigma_J)[1]
      }
    }
  } 
  return (max)
}

fwerfct <- function(m, Sigma, c, N, k){
  max = 0
  J <- indexfct(k, m)
  card_J <- length(J)
  Sigma_J <- subcovmatfct(J, card_J, Sigma)
  n <- 2^m-1
  if(card_J==1){
      max <- 1-pt(q=c, df=N-n)
  }
  else{
      max <- 1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1] 
  }
  return (max)
}

# Function for calculation of the PWER with the minimum number min
pwerfct_min <- function(m, piv, Sigma, c, N, indices, min){
  pwer = 0
  fac <- 1-length(indices)*min
  for(k in 1:(2^m-1)){
    J <- indexfct(k, m)
    card_J <- length(J)
    Sigma_J <- subcovmatfct(J, card_J, Sigma)
    n <- 2^m-1
    if(card_J==1){
      if(k %in% indices){
        pwer <- pwer + (min * (1-pt(q=c, df=N-n)))
      }
      else{
        pwer <- pwer + (fac*piv[k]* (1-pt(q=c, df=N-n)))
      }    
    }
    else{
      if(k %in% indices){
        pwer <- pwer + (min * (1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1])) 
      }
      else{
        pwer <- pwer + (fac*piv[k]*(1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1])) 
      }
    }
  } 
  return (pwer)
}

# Returns the index set J belonging to the k-th summand
indexfct <- function(k, m){
  v <- binvector(k, m)
  set <- vector()
  for(i in 1:(length(v))) {
    if(v[i]==1) set <- append(set, i)
  }  
  return (set)
}

# Conversion of the natural number k into a binary vector with d digits
binvector = function(k, d) {  
  vec = rev(as.numeric(intToBits(k)))
  return (vec[-(1:(length(vec)-d))])
}

# Returns the scale matrix
covmatfct <- function(piv, m){
  Sigma <- matrix(rep(1, m*m),nr=m)
  for(i in 1:m){
    for(j in 1:m){
      if(i!=j) Sigma[i,j] <- covfct(piv,i,j,m)
    }
  } 
  return(Sigma)
}

# Returns the scale matrix in case of a single tested treatment
covmatfct_s <- function(piv, m){
  Sigma <- matrix(rep(1, m*m),nr=m)
  for(i in 1:m){
    for(j in 1:m){
      if(i!=j) Sigma[i,j] <- covfct_s(piv,i,j,m)
    }
  } 
  return(Sigma)
}

# Returns the scale matrix in case of random treatment allocation numbers
covmatfct_ra <- function(m, S){
  Sigma <- matrix(rep(1, m*m),nr=m)
  for(i in 1:m){
    for(j in 1:m){
      if(i!=j) Sigma[i,j] <- covfct_ra(i,j,m,S)
    }
  } 
  return(Sigma)
}

# Creates the sub-correlation matrices
subcovmatfct <- function(J, card_J, Sigma){
  Sigma_J <- matrix(rep(1, card_J*card_J),nr=card_J)
  rho_J <- vector()
  for(i in J){
    for(j in J){
      if(i!=j) rho_J <- append(rho_J, Sigma[i,j])
    }
  }
  counter = 1
  for(i in 1:card_J){
    for(j in 1:card_J){
      if(i!=j) {
        Sigma_J[i,j] <- rho_J[counter] 
        counter <- counter+1
      }
    }
  }
  return(Sigma_J)
}

# Calculation of Cov(Z_i, Z_j) by formula (3)
covfct<- function(piv, i, j, m){
  pi_i <- 0
  pi_j <- 0
  U_i <- vector()
  U_j <- vector()
  for(k in 1:(2^m-1)){
    J <- indexfct(k,m)
    for(a in J){
      if(i==a) {
        pi_i <- pi_i + piv[k]*(length(J)+1)
        U_i <- append(U_i, k)
      }
      if(j==a){
        pi_j <- pi_j + piv[k]*(length(J)+1)
        U_j <- append(U_j, k)
      }
    }
  }
  pi_i_j <- 0
  for(b in U_i){
    for(c in U_j){
      if(b==c) pi_i_j <- pi_i_j + piv[b]*(length(indexfct(b,m))+1)
    }
  }
  return (ifelse(pi_i==0 | pi_j==0, 0, (pi_i_j)/(2*sqrt(pi_i * pi_j))))
}

# Calculation of Cov(Z_i, Z_j) in case of a single tested treatment
covfct_s<- function(piv, i, j, m){
  pi_i <- 0
  pi_j <- 0
  U_i <- vector()
  U_j <- vector()
  for(d in 1:(2^m-1)){
    for(e in indexfct(d,m)){
      if(i==e) {
        pi_i <- pi_i + piv[d] 
        U_i <- append(U_i, d)
      }
      if(j==e){ 
        pi_j <- pi_j + piv[d]
        U_j <- append(U_j, d)
      }
    }
  }
  pi_i_j <- 0
  for(f in U_i){
    for(g in U_j){
      if(f==g) pi_i_j <- pi_i_j + piv[f]
    }
  }
  return (ifelse(pi_i == 0 | pi_j == 0, 0, (pi_i_j)/sqrt(pi_i * pi_j)))
}

# Calculation of Cov(Z_i, Z_j) in case of random treatment allocation numbers
covfct_ra <- function(i, j, m, S){
  #First compute n_{i,C}, n_{j,C}, n_{i,T_i}, n_{j,T_j}
  n_iC= 0; n_jC=0; n_iTi= 0; n_jTj=0
  for(d in 1:(2^m-1)) if(i %in% indexfct(d,m)){
      n_iC <- n_iC + S[d,m+1]
      n_iTi <- n_iTi + S[d,i]
  }
  for(d in 1:(2^m-1)) if(j %in% indexfct(d,m)){
    n_jC <- n_jC + S[d,m+1]
    n_jTj <- n_jTj + S[d,j]
  }
  if(n_iC== 0 | n_jC==0| n_iTi== 0| n_jTj==0) return(0)
  #Numerator of formula from end of appendix B
  z <- 0
  for(d in 1:(2^m-1)) if(i %in% indexfct(d,m) & j %in% indexfct(d,m)){
    z <- z + S[d,m+1]/(n_iC * n_jC)
  }
  #Denominator 
  n1=0; n2=0
  for(d in 1:(2^m-1)) if(i %in% indexfct(d,m)) n1 <- n1 + S[d,i]/(n_iTi^2) + S[d,m+1]/(n_iC^2)
  for(d in 1:(2^m-1)) if(j %in% indexfct(d,m)) n2 <- n2 + S[d,j]/(n_jTj^2) + S[d,m+1]/(n_jC^2)
  return(z/sqrt(n1*n2))
}

# Calculation of the critical value c
critpwer <- function(m, piv, Sigma, alpha, N){
  res <- function(crit){pwerfct(m, piv, Sigma, crit, N)-alpha}
  return(uniroot(res, interval = c(1,5))$root)
}

# Calculation of the critical value c
critpwer_kv <- function(m, piv, Sigma, alpha){
  res <- function(crit){pwerfct_kv(m, piv, Sigma, crit)-alpha}
  return(uniroot(res, interval = c(1,5))$root)
}

# Calculation of the critical value c
critpwer_min <- function(m, piv, Sigma, alpha, N, indices, min){
  res <- function(crit){pwerfct_min(m, piv, Sigma, crit, N, indices, min)-alpha}
  return(uniroot(res, interval = c(1,5))$root)
}

# Calculates the relative prevalences from the biomarker probabilities
compute_pi <- function(m, p){
  pi <- rep(1, 2^m-1)
  for(k in 1:(2^m-1)){
    bin <- binvector(k,m)
    for(i in 1:m){
      if(bin[i]==1){
        pi[k] <- pi[k]*p[i]
      }
      else{
        pi[k] <- pi[k]*(1-p[i])
      }
    }
  }
  return (pi/sum(pi))
}

#Computes the vector of marginal sum estimators
marginal_pi <- function(m,piv){
  #Compute the vector of population prevalences
  pop_pi <- vector()
  for(i in 1:m){
    pi_i <- 0
    for(k in 1:(2^m-1)){
      J <- indexfct(k,m)
      if(i %in% J) pi_i <- pi_i + piv[k]
    }
    pop_pi <- append(pop_pi, pi_i)
  }
  #Compute two products for all J
  pi.est <- vector()
  for(k in 1:(2^m-1)){
    prod_1 <- 1
    prod_2 <- 1
    J <- indexfct(k,m)
    for(j in J) prod_1 <- prod_1 * pop_pi[j]
    for(j in 1:m) if(!(j %in% J)) prod_2 <- prod_2 * (1-pop_pi[j])
    pi.est <- append(pi.est, prod_1*prod_2)
  }
  return(pi.est)
}

#Matrix of sample sizes. For each row (representing the stratum given by the 
#indexfct), it defines allocation numbers under some randomly generated 
#allocation probabilities. 
allocation <- function(m,n){
  nrow <- 2^m-1
  ncol <- m+1
  
  # First generate random allocation probabilities
  aprobs <- matrix(0, nrow=nrow, ncol=ncol)
  for (i in 1:nrow) {
    tmts <- c(indexfct(i, m),m+1)
    aprobs[i, tmts] <- runif(length(tmts))
    aprobs[i, tmts] <- aprobs[i, tmts] / sum(aprobs[i, tmts])
  }
  
  #Equal allocation probabilities
  #aprobs <- matrix(0, nrow=nrow, ncol=ncol)
  #for (i in 1:nrow) {
  #  tmts <- c(indexfct(i, m),m+1)
  #  aprobs[i, tmts] <- 1/length(tmts)
  #}
  
  # Create the matrix with allocated values
  S <- matrix(0, nrow = nrow, ncol = ncol)
  
  for (i in 1:nrow) {
    tmts <- c(indexfct(i, m),m+1)
    # Allocate values based on allocation probabilities
    avalues <- rmultinom(n=1,size=n[i],prob=aprobs[i, tmts])
    # Fill the result matrix
    S[i,tmts] <- avalues
  }
  return(S)
}
