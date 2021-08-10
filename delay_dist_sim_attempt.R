###  Simulate delay_distributions
delay_dist_sim_updated<- function(df1){
  n.sim <- nrow(df1)
  expose.dist= rgamma(n.sim,shape=4,scale=1) #duration latent: Davies et al., The Lancet, 2020
  infect.dist= rgamma(n.sim,shape=4,scale=5/4) #duration of subclinical infectiousness -> temptative, this should be duration of clinical infectiousness, Davies, Lancet, 2020
  delay.dist= rgamma(n.sim,shape=2.34,scale=1/(2.59)) #duration infectiousness to test  #from Li.et al., Science  May 2020, estimates after 23rd Jan, parameters for Gamma should be a and b
  
  #First day is common across HHs and it is set as 15th of June
  #first.day.hh <-  startdate
  df1$date.exposed <- as.Date(df1$PCR_DATE - delay.dist - expose.dist,"%d-%m-%y")
  df1$date.infected.end <- as.Date(df1$PCR_DATE + infect.dist - delay.dist,"%d-%m-%y")
  df1$date.onset <- as.Date(df1$PCR_DATE - delay.dist,"%d-%m-%y")
  
  df1$first.date.hh <- min(df1$date.exposed[df1$infected==1], na.rm=T) #- 1
  df1$date.exposed[df1$infected==0] <- df1$first.date.hh[df1$infected==0]
  # Associate cal date to variables, equal to start
  NumDaysExp <- round(difftime(df1$date.exposed,startdate,units="days"))
  df1$day.exposed  <- as.numeric(NumDaysExp)     
  # # 
  NumDaysOnset <- round(difftime(df1$date.onset,startdate,units="days"))
  df1$day.infectious  <- as.numeric(NumDaysOnset)     
  # # 
  NumDaysInfEnd <- round(difftime(df1$date.infected.end,startdate,units="days"))
  df1$day.infectious.end   <- as.numeric(NumDaysInfEnd)     
  # 
  df1$day.infectious[df1$infected==0] <- NA
  df1$day.infectious.end[df1$infected==0] <- NA
  
  #df1$cal_date <- 1
  df1$first.day.hh <- min(df1$day.exposed, na.rm=T)
  df1$max.time.hh <- max(df1$day.infectious.end, na.rm=T)
  df1$enddate <- max(df1$date.infected.end, na.rm=T)+1
  return(df1)
}






