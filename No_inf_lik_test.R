##############
#Seed
##############
set.seed(8652)
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
library(lubridate)

source('./simulate_data.R')
source('./delay_dist_sim_attempt.R')
source('./data_manipulation_ksm_dose2.R')

##################
#Packages
##################

####################################################################
#Reading in Data
####################################################################

#For Josh:
#data_sim<-readRDS("C:/Users/jlw98/Desktop/df1.rds")

##### Load the data 
data_sim <- readRDS("./Data/simulated_data.rds")
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
data_spl <- lapply(sim.data.df.spl,function(x) x=suppressWarnings(delay_dist_sim_updated(x)))
df1 <- bind_rows(data_spl)
#######################################################
#Fake Data
#######################################################
first_day<-as.Date("2020-09-15")  #Just a guess
last_day<-as.Date("2021-04-20")  #Just a guess
cases<-rnorm(n = as.numeric(last_day - first_day + 1))


complete_dates<-rep(first_day,
                    times = (last_day - first_day + 1))
for(j in 1:(last_day - first_day + 1)){
  complete_dates[j]<-as.Date(first_day + j - 1)
}


#################################################################################################################################
#Likelihood Combinations and Calculations
#################################################################################################################################
params<-rnorm(5)  #Parameter values

df1$vax2dose_date[is.na(df1$vax2dose_date) == 1]<-last_day + 1  #Unvaccinated people (day after study ends)
df1$vax2dose_date<-floor_date(df1$vax2dose_date)
vax_cats<-unique(df1$vax2dose_date)
n_vax_cats<-length(vax_cats)

first_day<-floor_date(first_day)
last_day<-floor_date(last_day)
#last_day_HH<-floor_date(df1$first.date.hh)-1 #day before first exposure in the HH
complete_dates<-floor_date(complete_dates)

age1<-as.numeric(df1$AGE >= 10 & df1$AGE < 60)  #CHECK AGE GROUPS
age2<-as.numeric(df1$AGE >= 60)                            #CHECK AGE GROUPS

###############################################################################################################################################################
#Getting Counts:
#This should be calculated outside of the likelihood function (data preparation stage) and "s0, s1, s2, s0_sum, s1_sum, s2_sum" should be input to the function
###############################################################################################################################################################
s0<-matrix(0.00,
           nrow = n_vax_cats,
           ncol = (last_day - first_day + 1))
s1<-matrix(0.00,
           nrow = n_vax_cats,
           ncol = (last_day - first_day  + 1))
s2<-matrix(0.00,
           nrow = n_vax_cats,
           ncol = (last_day - first_day+ 1))
s0_sum<-matrix(0.00,
               nrow = n_vax_cats,
               ncol = (last_day - first_day+ 1))
s1_sum<-matrix(0.00,
               nrow = n_vax_cats,
               ncol = (last_day - first_day + 1))
s2_sum<-matrix(0.00,
               nrow = n_vax_cats,
               ncol = (last_day - first_day + 1))
for(j in 1:n_vax_cats){
  for(k in 1:(enddate - start + 1)){
    s0[j,k]<-sum((df1$vax2dose_date == vax_cats[j]) & (age1 == 0) & (age2 == 0) & ((df1$first.date.hh-1) == complete_dates[k])) ### check day before first exposure in the HH 
    s1[j,k]<-sum((df1$vax2dose_date == vax_cats[j]) & (age1 == 1) & ((df1$first.date.hh-1) == complete_dates[k]))                
    s2[j,k]<-sum((df1$vax2dose_date == vax_cats[j]) & (age2 == 1) & ((df1$first.date.hh-1) == complete_dates[k]))                
    
  }
}

for(j in 1:n_vax_cats){
  for(k in 1:(last_day- first_day)){
    s0_sum[j,k]<-sum(s0[j, (k + 1):(last_day - first_day + 1)])
    s1_sum[j,k]<-sum(s1[j, (k + 1):(last_day - first_day + 1)])
    s2_sum[j,k]<-sum(s2[j, (k + 1):(last_day - first_day + 1)])
    }
}

###############################################################################################################################################################
###############################################################################################################################################################

age0_log_probs<-rep(0.00,
                    times = n_vax_cats)
age1_log_probs<-rep(0.00,
                    times = n_vax_cats)
age2_log_probs<-rep(0.00,
                    times = n_vax_cats)
for(j in 1:n_vax_cats){
  #age0
  x<-cbind(1,  
           c(rep(0, times = as.numeric(vax_cats[j] - first_day)), rep(1, times = as.numeric(last_day - vax_cats[j] + 1))),
           0,
           0,
           cases)
  p0<-1.00/(1.00 + exp(-x%*%params))
  age0_log_probs[j]<-sum(s0_sum[j,]*log(1.00 - p0))
  
  #age1
  x<-cbind(1,
           c(rep(0, times = as.numeric(vax_cats[j] - first_day)), rep(1, times = as.numeric(last_day - vax_cats[j] + 1))),
           1,
           0,
           cases)
  p0<-1.00/(1.00 + exp(-x%*%params))
  age1_log_probs[j]<-sum(s1_sum[j,]*log(1.00 - p0))
  
  #age2
  x<-cbind(1,
           c(rep(0, times = as.numeric(vax_cats[j] - first_day)), rep(1, times = as.numeric(last_day - vax_cats[j] + 1))),
           0,
           1,
           cases)
  p0<-1.00/(1.00 + exp(-x%*%params))
  age2_log_probs[j]<-sum(s2_sum[j,]*log(1.00 - p0))
  
}
ll_contribution<-sum(age0_log_probs +
                       age1_log_probs +
                       age2_log_probs)

















