
library(multcomp)
library(randcorr)

# Function for calculation of the PWER
#m: number of populations
#pi: prevalence vector
#ALLOC: sample size matrix
#c: critical value
#distribution: "t" for t-distr, "z" for z-distr 
#min: minimum prevalence for empty strata
pwer <- function(m, pi, ALLOC, c, distr="t", min=0){
  pwer = 0
  n <- 2^m-1 # number of strata
  N <- round(sum(ALLOC))
  Sigma <- covmat(m, ALLOC)
  pi.est <- apply(ALLOC, 1, sum)/N
  indices <- which(pi.est < min) #check which strata are smaller than min
  if(length(indices)==0) fac <- 1 
  else fac <- (1-length(indices)*min)/(1-sum(pi.est[indices]))
  for(k in 1:n){
    J <- indexfct(k, m)
    card_J <- length(J)
    Sigma_J <- subcovmat(J, card_J, Sigma)
    if(card_J==1 & distr=="t" & k %in% indices){
      pwer <- pwer + (min * (1-pt(q=c, df=N-n)))
    }
    if(card_J==1 & distr=="t" & !(k %in% indices)){
      pwer <- pwer + (pi[k]* fac * (1-pt(q=c, df=N-n)))
    }
    if(card_J==1 & distr=="z" & k %in% indices){
      pwer <- pwer + (min * (1-pnorm(c)))
    }
    if(card_J==1 & distr=="z" & !(k %in% indices)){
      pwer <- pwer + (pi[k]* fac * (1-pnorm(c)))
    }
    if(card_J>1 & distr=="t" & k %in% indices){
      pwer <- pwer + (min* (1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1])) 
    }
    if(card_J>1 & distr=="t" & !(k %in% indices)){
      pwer <- pwer + (pi[k] * fac * (1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1])) 
    }
    if(card_J>1 & distr=="z" & k %in% indices){
      pwer <- pwer + (min * (1-pmvnorm(upper=rep(c,card_J), corr=Sigma_J)[1]))  
    }
    if(card_J>1 & distr=="z" & !(k %in% indices)){
      pwer <- pwer + (pi[k] * fac * (1-pmvnorm(upper=rep(c,card_J), corr=Sigma_J)[1]))  
    }
  }
  return (pwer)
}


# Function for calculation of the maximal strata-wise FWER
max_fwer <- function(m, ALLOC, c, distr="t"){
  max = 0
  n <- 2^m-1 
  N <- round(sum(ALLOC))
  Sigma <- covmat(m, ALLOC)
  for(k in 1:n){
    J <- indexfct(k, m)
    card_J <- length(J)
    Sigma_J <- subcovmat(J, card_J, Sigma)
    if(card_J==1 & distr=="t"){
      if(1-pt(q=c, df=N-n)>max) max <- 1-pt(q=c, df=N-n)
    }
    if(card_J==1 & distr=="z"){
      if(1-pnorm(c)>max) max <- 1-pnorm(c)
    }
    if(card_J>1 & distr=="t"){
      if(1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1]>max){
        max <- 1-pmvt(upper=rep(c,card_J), sigma=Sigma_J, df=N-n)[1] 
      }
    }
    if(card_J>1 & distr=="z"){
      if(1-pmvnorm(upper=rep(c,card_J), corr=Sigma_J)[1]>max){
        max <- 1-pmvnorm(upper=rep(c,card_J), corr=Sigma_J)[1]
      }
    }
  } 
  return (max)
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

# Computes the scale matrix
covmat <- function(m, ALLOC){
  Sigma <- matrix(rep(1,m*m),nr=m)
  for(i in 1:m){
    for(j in 1:m){
      if(i!=j) Sigma[i,j] <- cov(i,j,m,ALLOC)
    }
  } 
  return(Sigma)
}


# Calculation of Cov(Z_i, Z_j) for i neq j
# ALLOC: matrix of allocation numbers
cov <- function(i, j, m, ALLOC){ 
  if(ncol(ALLOC)!=2){ #case of pairwise different treatments
    #First compute n_{i,C}, n_{j,C}, n_{i,T_i}, n_{j,T_j}
    n_iC=0; n_jC=0; n_iTi=0; n_jTj=0
    for(d in 1:(2^m-1)) if(i %in% indexfct(d,m)){
      n_iC <- n_iC + ALLOC[d,m+1]
      n_iTi <- n_iTi + ALLOC[d,i]
    }
    for(d in 1:(2^m-1)) if(j %in% indexfct(d,m)){
      n_jC <- n_jC + ALLOC[d,m+1]
      n_jTj <- n_jTj + ALLOC[d,j]
    }
    if(n_iC== 0 | n_jC==0| n_iTi== 0| n_jTj==0) return(0)
    #Compute V_i and V_j
    V_i=0; V_j=0
    for(d in 1:(2^m-1)) if(i %in% indexfct(d,m)) V_i <- V_i + ALLOC[d,i]/(n_iTi^2) + ALLOC[d,m+1]/(n_iC^2)
    for(d in 1:(2^m-1)) if(j %in% indexfct(d,m)) V_j <- V_j + ALLOC[d,j]/(n_jTj^2) + ALLOC[d,m+1]/(n_jC^2)
    
    #Numerator of formula from end of appendix B
    retv <- 0
    for(d in 1:(2^m-1)) if(i %in% indexfct(d,m) & j %in% indexfct(d,m)){
      retv <- retv + ALLOC[d,m+1]/(n_iC * n_jC)
    }
    return(retv/sqrt(V_i*V_j))
  } 
  else{ # single treatment case
    n_iC=0; n_jC=0; n_iT=0; n_jT=0
    for(d in 1:(2^m-1)) if(i %in% indexfct(d,m)){
      n_iC <- n_iC + ALLOC[d,2]
      n_iT <- n_iT + ALLOC[d,1]
    }
    for(d in 1:(2^m-1)) if(j %in% indexfct(d,m)){
      n_jC <- n_jC + ALLOC[d,2]
      n_jT <- n_jT + ALLOC[d,1]
    }
    if(n_iC== 0 | n_jC==0| n_iT== 0| n_jT==0) return(0)
    V_i=0; V_j=0
    for(d in 1:(2^m-1)) if(i %in% indexfct(d,m)) V_i <- V_i + ALLOC[d,1]/(n_iT^2) + ALLOC[d,2]/(n_iC^2)
    for(d in 1:(2^m-1)) if(j %in% indexfct(d,m)) V_j <- V_j + ALLOC[d,1]/(n_jT^2) + ALLOC[d,2]/(n_jC^2)
    retv <- 0
    for(d in 1:(2^m-1)) if(i %in% indexfct(d,m) & j %in% indexfct(d,m)){
      retv <- retv + ALLOC[d,1]/(n_iT*n_jT) + ALLOC[d,2]/(n_iC * n_jC)
    }
    return(retv/sqrt(V_i*V_j))
  }
}


# Creates the sub-correlation matrices
subcovmat <- function(J, card_J, Sigma){
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


# Calculation of the critical value c
critpwer <- function(m, pi, ALLOC, alpha, distr="t", min="0"){
  res <- function(c){pwer(m, pi, ALLOC, c, distr, min)-alpha}
  return(uniroot(res, interval = c(0,20))$root)
}


# Calculates the prevalences of the relevant strata
compute_pi <- function(m, tau){
  tau[2:(2^m)]/(1-tau[1])
}

# Calculates the prevalences (including the stratum where none of the biomarkers
# are present)
compute_tau <- function(m, p){
  tau <- rep(1,2^m)
  for(k in 1:(2^m)){
    bin <- binvector(k-1,m)
    for(i in 1:m){
      if(bin[i]==1){
        tau[k] <- tau[k]*p[i]
      }
      else{
        tau[k] <- tau[k]*(1-p[i])
      }
    }
  }
  return(tau)
}


#Computes the vector of marginal sum estimators
marginal_pi <- function(m,tau){
  #Compute the marginal prevalences
  p.est <- rep(0,m)
  for(i in 1:m){
    for(k in 1:(2^m-1)){
      J <- indexfct(k,m)
      if(i %in% J) p.est[i] <- p.est[i] + tau[k+1]
    }
  }
  tau.marg <- compute_tau(m,p.est)
  compute_pi(m,tau.marg)
}

#Creates a matrix containing the allocation numbers. The rows correspond to the 
#strata (given by the indexfct) and the columns to the treatments (control in 
#last column)
#n: vector of strata-wise sample sizes
#case: "s" for single treatment case, "d" for pairwise different treatments
#alloc: "equal" for equal allocation probabilities within the strata, 
#       "random" for random allocation numbers based on equal allocation probabilities.
#       "randomp" for random allocation numbers based on random allocation probabilities.
allocation <- function(m,n,alloc="equal",case="d"){
  nrow <- 2^m-1
  if(case=="d") ncol <- m+1
  if(case=="s") ncol <- 2
  ALLOC <- matrix(0, nrow = nrow, ncol = ncol) #matrix with allocation numbers
  if(alloc=="equal"){
    for (i in 1:nrow) {
      if(case=="d") tmts <- c(indexfct(i,m),m+1)
      if(case=="s") tmts <- c(1,2)
      ALLOC[i, tmts] <- n[i]*1/length(tmts)
    }
    return(ALLOC)
  }
  if(alloc=="random"){
    aprobs <- matrix(0, nrow=nrow, ncol=ncol) #matrix with allocation probabilities
    for (i in 1:nrow) {
      if(case=="d") tmts <- c(indexfct(i,m),m+1)
      if(case=="s") tmts <- c(1,2)
      aprobs[i, tmts] <- 1/length(tmts)
    }
  }
  if(alloc=="randomp"){
    aprobs <- matrix(0, nrow=nrow, ncol=ncol)
    for (i in 1:nrow) {
      if(case=="d") tmts <- c(indexfct(i,m),m+1)
      if(case=="s") tmts <- c(1,2)
      aprobs[i, tmts] <- runif(length(tmts))
      aprobs[i, tmts] <- aprobs[i, tmts] / sum(aprobs[i, tmts])
    }
  }
  # Generate allocation numbers based on allocation probabilities
  for (i in 1:nrow) {
    if(case=="d") tmts <- c(indexfct(i,m),m+1)
    if(case=="s") tmts <- c(1,2)
    avalues <- rmultinom(n=1,size=n[i],prob=aprobs[i, tmts])
    ALLOC[i,tmts] <- avalues
  }
  return(ALLOC)
}
