# Check timings to evaluate log-likelihood
library(data.table)
source('./Chain_bin_lik.R')
# Load 1 dataset among the different iterations of the delay distributions
ds <- readRDS(file=paste("./LatentData_primarymodel/latentData_",1,".rds"))


# Set values for the parameters
alpha0 = -3.94
delta0 = runif(1,-8.25,-6)
beta2 = -0.9
kappa2 = -0.7
gamma1 = 0.6
gamma2 = 1
delta1 = runif(1,0.6,0.9)
params <- c(alpha0,delta0,beta2,kappa2,gamma1,gamma2,delta1)

ptm <- proc.time()

chain_bin_lik(params,ds$Y,ds$X)

proc.time() - ptm

