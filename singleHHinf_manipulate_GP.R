######################################################
#Fake Data
######################################################
singleHHinf_manipulate <- function(df){

  start<-floor_date(startdate)
  enddate<- floor_date(enddate)
  complete_dates<-start:enddate
  
  #################################################################################################################################
  #Likelihood Combinations and Calculations
  #################################################################################################################################
  df$FIRST_VACCINE_DATE <- ifelse((df$FIRST_VACCINE_DATE+7)<=as.character(enddate),as.character(df$FIRST_VACCINE_DATE+7),as.character(enddate+1))
  df$FIRST_VACCINE_DATE[is.na(df$FIRST_VACCINE_DATE) == 1]<-as.character(enddate + 1)  #Unvaccinated people (day after study ends)
  df$FIRST_VACCINE_DATE<-floor_date(as.Date(df$FIRST_VACCINE_DATE))
  
  df$SECOND_VACCINE_DATE <- ifelse((df$SECOND_VACCINE_DATE+10)<=as.character(enddate),as.character(df$SECOND_VACCINE_DATE+10),as.character(enddate+1))
  df$SECOND_VACCINE_DATE[is.na(df$SECOND_VACCINE_DATE) == 1]<-as.character(enddate + 1)  #Unvaccinated people (day after study ends)
  df$SECOND_VACCINE_DATE<-floor_date(as.Date(df$SECOND_VACCINE_DATE))
    
  #vax_cats<-unique(df$FIRST_VACCINE_DATE) # I think it might mess things up if there are days after vaccination started in which nobody is reported to be vaccinated
  #n_vax_cats<-length(vax_cats)
    
  df$exposed.date <- floor_date(df$PCR_DATE - rgamma(nrow(df),shape=4,scale=1)-rgamma(nrow(df), shape=2.34,scale=1/(2.59)), unit='day') #note: need to use floor date otherwise have fractional days, which causes problems later
    
  vax1 <- matrix(0.00, nrow = length(df), ncol = (enddate - start + 1))
  vax2 <- matrix(0.00, nrow = length(df), ncol = (enddate - start + 1))
  vax2_wane <- matrix(0.00, nrow = length(df), ncol = (enddate - start + 1))
  suscept <- matrix(0.00, nrow = length(df), ncol = (enddate - start + 1))
  infect <- matrix(0.00, nrow = length(df), ncol = (enddate - start + 1))
  for (j in 1:(enddate - start + 1)){
    vax1[,j] <- as.numeric((df$FIRST_VACCINE_DATE - start + 1) < j)
    vax2[,j] <- as.numeric((df$SECOND_VACCINE_DATE - start + 1) < j)
    vax2_wane[,j] <- as.numeric((df$SECOND_VACCINE_DATE - start + 121) < j)
    suscept[,j] <- as.numeric((df$exposed.date - start + 1) > j)
    infect[,j] <- as.numeric((df$exposed.date - start + 1) == j)
  }
  #vax2 <- vax[,(vax_cats-start)] # only need to keep vax status for unique vax dates
  #vax2_wane <- vax2_wane[,(vax_cats-start)] # only need to keep vax status for unique vax dates
  vax1_part <- vax1-vax2 # partially vaccinated individuals
  vax2_new <- vax2-vax2_wane # individuals who were vaccinated <180 days ago
  
  # Using logical indexing below
  age0<-(df$AGE < 12)                 #CHECK AGE GROUPS
  age1<-(df$AGE >= 12 & df$AGE < 60)  #CHECK AGE GROUPS
  age2<-(df$AGE >= 60)                #CHECK AGE GROUPS
  
  month<-matrix(rep(0,13*(enddate - start + 1)), nrow = (enddate - start + 1), ncol = 14)
  month[1:16,1]<-rep(1,16) #Jun 15-30, 2020
  month[17:47,2]<-rep(1,31) #Jul 1-31
  month[48:78,3]<-rep(1,31) #Aug 1-31
  month[79:108,4]<-rep(1,30) #Sep 1-30
  month[109:139,5]<-rep(1,31) #Oct 1-31
  month[140:169,6]<-rep(1,30) #Nov 1-30
  month[170:200,7]<-rep(1,31) #Dec 1-31
  month[201:231,8]<-rep(1,31) #Jan 1-31, 2021
  month[232:259,9]<-rep(1,28) #Feb 1-28
  month[260:290,10]<-rep(1,31) #Mar 1-31
  month[291:320,11]<-rep(1,30) #Apr 1-30
  month[321:351,12]<-rep(1,31) #May 1-31
  month[352:381,13]<-rep(1,30) #June 1-30
  month[382:409,14]<-rep(1,28) #July 1-28
  
  dv<-c(rep(0,351), rep(1,58)) #Delta variant period
    
  ###############################################################################################################################################################
  #Getting Counts:
  #This should be calculated outside of the likelihood function (data preparation stage) 
  ###############################################################################################################################################################
  
  # Number vax on day i who were infected on day j
  i0vax <- t(vax2_new[age0,])%*%infect[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  i1vax <- t(vax2_new[age1,])%*%infect[age1,]
  i2vax <- t(vax2_new[age2,])%*%infect[age2,]
  
  i0wane <- t(vax2_wane[age0,])%*%infect[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  i1wane <- t(vax2_wane[age1,])%*%infect[age1,]
  i2wane <- t(vax2_wane[age2,])%*%infect[age2,]
  
  i0pvax <- t(vax1_part[age0,])%*%infect[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  i1pvax <- t(vax1_part[age1,])%*%infect[age1,]
  i2pvax <- t(vax1_part[age2,])%*%infect[age2,]
  
  i0unvax <- t(1-vax1[age0,])%*%infect[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  i1unvax <- t(1-vax1[age1,])%*%infect[age1,]
  i2unvax <- t(1-vax1[age2,])%*%infect[age2,]
  
  # Number vax on day i who were susceptible on day j
  s0vax <- t(vax2_new[age0,])%*%suscept[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  s1vax <- t(vax2_new[age1,])%*%suscept[age1,]
  s2vax <- t(vax2_new[age2,])%*%suscept[age2,]
  
  s0wane <- t(vax2_wane[age0,])%*%suscept[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  s1wane <- t(vax2_wane[age1,])%*%suscept[age1,]
  s2wane <- t(vax2_wane[age2,])%*%suscept[age2,]
  
  s0pvax <- t(vax1_part[age0,])%*%suscept[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  s1pvax <- t(vax1_part[age1,])%*%suscept[age1,]
  s2pvax <- t(vax1_part[age2,])%*%suscept[age2,]
  
  s0unvax <- t(1-vax1[age0,])%*%suscept[age0,] # this should be matrix multiplication of (n_vax_cat x n_indiv) x (n_indiv x tot_days) = (n_vax_cat x tot_days) 
  s1unvax <- t(1-vax1[age1,])%*%suscept[age1,]
  s2unvax <- t(1-vax1[age2,])%*%suscept[age2,]
  
  for (j in 1:14){
    pt0vax[j] <- sum(s0vax%*%month[,j]) #Susceptible person-time among vaccinated individuals in month j
    pt1vax[j] <- sum(s1vax%*%month[,j])
    pt2vax[j] <- sum(s2vax%*%month[,j])
    
    pt0wane[j] <- sum(s0wane%*%month[,j]) #Susceptible person-time among vaccinated individuals (>120d) in month j
    pt1wane[j] <- sum(s1wane%*%month[,j])
    pt2wane[j] <- sum(s2wane%*%month[,j])
       
    pt0pvax[j] <- sum(s0pvax%*%month[,j]) #Susceptible person-time among partially vaccinated individuals in month j
    pt1pvax[j] <- sum(s1pvax%*%month[,j])
    pt2pvax[j] <- sum(s2pvax%*%month[,j]) 
    
    pt0unvax[j] <- sum(s0unvax%*%month[,j]) #Susceptible person-time among unvaccinated individuals in month j
    pt1unvax[j] <- sum(s1unvax%*%month[,j])
    pt2unvax[j] <- sum(s2unvax%*%month[,j])   
    
    infect0vax[j] <- sum(i0vax%*%month[,j]) #Number of vaccinated individuals infected in month j
    infect1vax[j] <- sum(i1vax%*%month[,j])
    infect2vax[j] <- sum(i2vax%*%month[,j])
    
    infect0wane[j] <- sum(i0wane%*%month[,j]) #Number of vaccinated individuals (>120d) infected in month j
    infect1wane[j] <- sum(i1wane%*%month[,j])
    infect2wane[j] <- sum(i2wane%*%month[,j])
   
    infect0pvax[j] <- sum(i0pvax%*%month[,j]) #Number of partially vaccinated individuals infected in month j
    infect1pvax[j] <- sum(i1pvax%*%month[,j])
    infect2pvax[j] <- sum(i2pvax%*%month[,j])
 
    infect0unvax[j] <- sum(i0unvax%*%month[,j]) #Number of unvaccinated individuals infected in month j
    infect1unvax[j] <- sum(i1unvax%*%month[,j])
    infect2unvax[j] <- sum(i2unvax%*%month[,j])
  }
  
  #out.list <- list('s0'=s0,'s1'=s1,'s2'=s2,'s0_sum'=s0_sum,'s1_sum'=s1_sum,'s2_sum'=s2_sum,'vax_cats'=vax_cats,'n_vax_cats'=n_vax_cats)
  out.list <- list('dv'=dv, 'month'=month, 'pt0vax'=pt0vax, 'pt1vax'=pt1vax, 'pt2vax'=pt2vax, 'pt0unvax'=pt0unvax, 'pt1unvax'=pt1unvax, 'pt2unvax'=pt2unvax, 
                   'pt0wane'=pt0wane, 'pt1wane'=pt1wane, 'pt2wane'=pt2wane, 'pt0pvax'=pt0pvax, 'pt1pvax'=pt1pvax, 'pt2pvax'=pt2pvax,  
                   'infect0vax'=infect0vax, 'infect1vax'=infect1vax, 'infect2vax'=infect2vax, 'infect0unvax'=infect0unvax, 'infect1unvax'=infect1unvax, 'infect2unvax'=infect2unvax,
                   'infect0wane'=infect0wane, 'infect1wane'=infect1wane, 'infect2wane'=infect2wane, 'infect0pvax'=infect0pvax, 'infect1pvax'=infect1pvax, 'infect2pvax'=infect2pvax)
                   
  return(out.list)
}


