######################################################
#Fake Data
######################################################
singleHHinf_manipulate_GP <- function(df){

  start<-floor_date(startdate)
  enddate<- floor_date(enddate)
  complete_dates<-start:enddate
  
  #################################################################################################################################
  #Likelihood Combinations and Calculations
  #################################################################################################################################
  df$SECOND_VACCINE_DATE <- ifelse((df$SECOND_VACCINE_DATE+10)<=as.character(enddate),as.character(df$SECOND_VACCINE_DATE+10),as.character(enddate+1))
  df$SECOND_VACCINE_DATE[is.na(df$SECOND_VACCINE_DATE) == 1]<-as.character(enddate + 1)  #Unvaccinated people (day after study ends)
  df$SECOND_VACCINE_DATE<-floor_date(as.Date(df$SECOND_VACCINE_DATE))
  vax_cats<-unique(df$SECOND_VACCINE_DATE)
  n_vax_cats<-length(vax_cats)
    
  df$exposed.date <- floor_date(df$PCR_DATE - rgamma(nrow(df),shape=4,scale=1)-rgamma(nrow(df), shape=2.34,scale=1/(2.59)), unit='day') #note: need to use floor date otherwise have fractional days, which causes problems later
    
  vax <- matrix(0.00, nrow = length(df), ncol = (enddate - start + 1))
  suscept <- matrix(0.00, nrow = length(df), ncol = (enddate - start + 1))
  infect <- matrix(0.00, nrow = length(df), ncol = (enddate - start + 1))
  for (j in 1:(enddate - start + 1)){
    vax[,j] <- as.numeric((df$SECOND_VACCINE_DATE - start + 1) < j)
    suscept[,j] <- as.numeric((df$exposed.date - start + 1) > j)
    infect[,j] <- as.numeric((df$exposed.date - start + 1) == j)
  }
  vax2 <- vax[,(vax_cats-start)] # only need to keep vax status for unique vax dates
  
  # Using logical indexing below
  age0<-(df$AGE < 12)                 #CHECK AGE GROUPS
  age1<-(df$AGE >= 12 & df$AGE < 60)  #CHECK AGE GROUPS
  age2<-(df$AGE >= 60)                #CHECK AGE GROUPS
  
  ###############################################################################################################################################################
  #Getting Counts:
  #This should be calculated outside of the likelihood function (data preparation stage) and "s0, s1, s2, s0_sum, s1_sum, s2_sum" should be input to the function
  ###############################################################################################################################################################
  
  # Number vax on day i who were infected on day j
  s0 <- t(vax2[age0,])%*%infect[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  s1 <- t(vax2[age1,])%*%infect[age1,]
  s2 <- t(vax2[age2,])%*%infect[age2,]
  
  # Number vax on day i who were susceptible on day j
  s0_sum <- t(vax2[age0,])%*%suscept[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  s1_sum <- t(vax2[age1,])%*%suscept[age1,]
  s2_sum <- t(vax2[age2,])%*%suscept[age2,]
  
  
  out.list <- list('s0'=s0,'s1'=s1,'s2'=s2,'s0_sum'=s0_sum,'s1_sum'=s1_sum,'s2_sum'=s2_sum,'vax_cats'=vax_cats,'n_vax_cats'=n_vax_cats)
  return(out.list)
}


