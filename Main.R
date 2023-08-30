
m <- 2:8
Nsim <- 10^4
alpha <- 0.025
N <- 500

source("Exportfunctions.R")

######################### True PWER #############################

source("Simulations.R")

#Create dataframes with simulated values in the different cases of section 4
#Execution for different m, with N=500 fixed
for(i in 1:8) {
  dfi <- paste0("df", i)
  df <- sapply(m, sim, Nsim=Nsim, alpha=alpha, N=N, case=i, value="pwer")
  df <- as.data.frame(df)
  names(df) <- paste("m=", as.character(m), sep="")
  assign(dfi, df)
}

# Boxplots
for(i in c(1,2,6,7,8)){dfi <- paste0("df", i); boxplots(get(dfi), i, c(0.024,0.026),"pwer")}
for(i in c(3,5)){dfi <- paste0("df", i); boxplots(get(dfi), i, c(0.02375,0.02625), "pwer")}
boxplots(df4, 4, c(0.02275,0.02725), "pwer")

# Summary tables
for(i in 1:8){dfi <- paste0("df", i); summary_table(get(dfi), i)}


# Execution for different N, with m=3,4 fixed
N <- c(25,50,100,150,200,500)
for(i in 9:10){
  dfi <- paste0("df", i)
  df <- sapply(N, sim, m=i-6, Nsim=Nsim, alpha=alpha, case=1, value="pwer")
  df <- as.data.frame(df)
  names(df) <- paste("N=", as.character(N), sep="")
  assign(dfi, df)
}

# Boxplots
for(i in 9:10){dfi <- paste0("df", i); boxplotsN(get(dfi), i, c(0.02,0.03),"pwer")}

#Summary tables
for(i in 9:10){dfi <- paste0("df", i); summary_table(get(dfi), i)}


####################### Maximal strata-wise FWER #########################
N <- 500
#Create dataframes with simulated values
for(i in 11:18) {
  dfi <- paste0("df", i)
  df <- sapply(m, sim, Nsim=Nsim, alpha=alpha, N=N, case=i-10, value="fwer")
  df <- as.data.frame(df)
  names(df) <- paste("m=", as.character(m), sep="")
  assign(dfi, df)
}

# Boxplots
for(i in 11:18){dfi <- paste0("df", i); boxplots(get(dfi), i, c(0.024,0.08),"fwer")}

#Summary tables
for(i in 11:18){dfi <- paste0("df", i); summary_table(get(dfi), i)}

# Execution for different N, with m=3,4 fixed
N <- c(25,50,100,150,200,500)
for(i in 19:20){
  dfi <- paste0("df", i)
  df <- sapply(N, sim, m=i-16, Nsim=Nsim, alpha=alpha, case=1, value="fwer")
  df <- as.data.frame(df)
  names(df) <- paste("N=", as.character(N), sep="")
  assign(dfi, df)
}

# Boxplots
boxplotsN(df19, 19, c(0.024,0.075),"fwer")
boxplotsN(df20, 20, c(0.024,0.08),"fwer")

#Summary tables
for(i in 19:20){dfi <- paste0("df", i); summary_table(get(dfi), i)}



######################## Minimal number ##############################
source("min.R")
N <- 500

for(i in 21:27) {
  dfi <- paste0("df", i)
  df <- sim_min(m=i-19, Nsim=Nsim, alpha=alpha, N=N)
  df <- as.data.frame(df)
  assign(dfi, df)
}

# Boxplots
for(i in 21:27){dfi <- paste0("df", i); boxplots_min(get(dfi), i)}

#Summary tables
for(i in 21:27){dfi <- paste0("df", i); summary_table_min(get(dfi), i, Nsim)}


