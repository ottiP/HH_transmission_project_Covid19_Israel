###### Call to the libraries
library(boot)
library(dplyr)
library(pbapply)
library(matrixStats)
library(reshape2)
library(data.table)
library(stats4)
library(foreach)
library(parallel)
library(msm)

source('./simulate_data.R')
source('./delay_dist_sim_updated.R')
source('./Chain_bin_lik.R')
source('./data_manipulation_ksm_dose2.R')


#Helper functions for processing and calling the data

gen.latent.data.par <- function(ds,i,cl1){
  latentdata <- pblapply(ds,delay.gen.ksm.singleHH.age,cl=cl1)
  X <- lapply(latentdata, function(x) x<- x$X)
  Y <- lapply(latentdata, function(x) x<- x$Y)
  X <- do.call('rbind.data.frame',X)
  Y <- unlist(Y)
  out.list <- list('X'=X,'Y'=Y)
  saveRDS(out.list,file=paste("../LatentData_primarymodel/latentData_",i,".rds"),compress=TRUE)
  X <- NULL
  Y <- NULL
  latentdata <- NULL
}

parRead <- function(i){
  ds <- readRDS(file=paste("../LatentData_primarymodel/latentData_",i,".rds"))
  output <- model.run(ds$X,ds$Y)
  ds <- NULL
  return(output)
}

model.run <- function(X,Y){
  ### Initiate parameter values with estimates derived from analysis on subset of HH and on HH with single occupants
  alpha0 = -3.94
  delta0 = runif(1,-8.25,-6)
  beta2 = -0.9
  kappa2 = -0.7
  gamma1 = 0.6
  gamma2 = 1
  delta1 = runif(1,0.6,0.9)
  params <- c(alpha0,delta0,beta2,kappa2,gamma1,gamma2,delta1)
  nlm(chain_bin_lik,p=params,Y=Y,X=X,hessian=TRUE)
}

##### Load the data 
data_sim <- readRDS("../Data/simulated_data.rds")
startdate <- as.Date("2020-08-29") ### allow for up to 17 days prior to the start of the study period to infections to occur
start <- as.Date("2020-09-15")
enddate <- as.Date("2021-03-24")#as.Date(max(data_sim$PCR_DATE,na.rm = TRUE),"%d-%m-%y")
pos_test_complete <- table(data_sim$PCR_DATE[data_sim$infected==1])
data_spl <- split(data_sim,data_sim$HH_CERTAIN)
data_spl_restricted <- lapply(data_spl,function(x) x=x[(min(x$PCR_DATE,na.rm = TRUE)>=start)&(max(x$PCR_DATE,na.rm=TRUE)<=enddate),])
data_restricted <- bind_rows(data_spl_restricted)
data_restricted <- data_restricted[rowSums(is.na(data_restricted))!=ncol(data_restricted),]
  
#### First part: manipulate the data 
sim.data.df.spl <- split(data_restricted,data_restricted$HH_CERTAIN)


fileNumber <- 1


nCores <- detectCores()-1
cl1 <- makeCluster(nCores)

clusterEvalQ(cl1, {
  library(lme4,quietly=TRUE)
  library(data.table,quietly=TRUE)
  library(dplyr,quietly=TRUE)
  library(mltools,quietly=TRUE)
  library(pbapply,quietly=TRUE)
})

clusterExport(cl1,c('sim.data.df.spl','delay_dist_sim_updated_age',
                    'gen.latent.data.par','start','startdate','enddate','fileNumber','pos_test_complete'),environment())
  
  for(i in 1:fileNumber){
    gen.latent.data.par(ds=sim.data.df.spl,i,cl1)
    if(i%%10==0){
      gc()
    }
  }

stopCluster(cl1)

#### Second part: fit the model using the manipulated data

cl1 <- makeCluster(nCores)
clusterEvalQ(cl1, {
  library(lme4,quietly=TRUE)
  library(data.table,quietly=TRUE)
  library(dplyr,quietly=TRUE)
  library(mltools,quietly=TRUE)
  library(pbapply,quietly=TRUE)
  require(foreach)
  require(doParallel)
  registerDoParallel(cores=4)
})

clusterExport(cl1,c('model.run','parRead','chain_bin_lik','fileNumber'),envir=environment())

out.list <- pblapply(seq(1,fileNumber),parRead,cl=cl1)
saveRDS(out.list,file="../Results_primarymodel/results_primary.rds")
stopCluster(cl1)


#### Analysis of the results





table((data_sim$PCR_DATE[data_sim$infected==1]))
