### preprocessing on the data: 
set.seed(1234)
library(reshape2)

###Simulate  people
gen.hh <- function(idN, CPI=(1-0.90), prob.trans.day=(1-0.90), prop.vax1=0.5, prop.vax2=0.75, irr.vax1=1,irr.vax2=1, IRR.comm=1){
  
  HH.size <- min(2+ rpois(n=1,1.5),8) #cap at 8; need at least 2 people
  df1 <- as.data.frame(matrix(NA, nrow=HH.size, ncol=2))
  names(df1) <- c('ID', 'HH_CERTAIN')
  df1$HH_CERTAIN <- idN
  df1$ID <- 1:HH.size
  
  df1$AGE <- round(runif(nrow(df1),min=0,max=85),0)
  
  df1$vax1dose <- rbinom(nrow(df1), 1,prop.vax1)
  df1$vax2dose <- rbinom(nrow(df1), 1,prop.vax2)*df1$vax1dose
  
  sample.date <- sample(seq(as.Date('2020-12-20'), as.Date('2021-03-24'), by = "day"), length(df1$vax1dose[df1$vax1dose==1]))
  df1$vax1dose_date <- NA
  df1$vax1dose_date[df1$vax1dose==1] <- as.character(sample.date)
  
  df1$vax2dose_date <- NA*df1$vax2dose
  df1$vax2dose_date[df1$vax2dose==1] <- as.character(as.Date(df1$vax1dose_date[df1$vax2dose==1]) + round(runif(length(df1$vax1dose[df1$vax2dose==1]),min=21,max=28),0))
  
  df1$vax1dose_date <- as.Date(df1$vax1dose_date)
  df1$vax2dose_date <- as.Date(df1$vax2dose_date)
  #These are time distributions that specify for the person if they are infected, how much time to add before infectious, and how long infectious
  expose.dist= rgamma(length(df1$ID),shape=4,scale=1) #duration latent
  infect.dist= rgamma(length(df1$ID),shape=4,scale=5/4) #duration infectiousness
 
   exposed.status <- matrix(NA, nrow=nrow(df1), ncol=200)
  infect.status <- matrix(NA, nrow=nrow(df1), ncol=200)
  n.infect.prev <- matrix(NA, nrow=nrow(df1), ncol=200)
  
  prob.infect.day <- prob.trans.day * irr.vax1^df1$vax1dose*irr.vax2^df1$vax2dose  #prob of being infected per day, per exposure,
  prob.uninfect.day <- 1 - prob.infect.day
  prob.uninf.day.comm <- 1 - CPI * irr.vax1^df1$vax1dose*irr.vax2^df1$vax2dose  
  
  #exposed.status[,1] <- df1$index*rbinom(nrow(df1),1, CPI*IRR.comm^(df1$vax1dose)) #10% chance that a tested index is positive 
  infect.status[,1] <- 0
  n.infect.prev[,1] <- 0
  exposed.status[,1] <- 1 - rbinom(nrow(df1), 1, prob.uninf.day.comm )  #exponent ensure once you are exposed, you stay in that category
  
  
  for( i in 2:ncol(infect.status)){
    day.exposed <- apply(exposed.status,1, function(x) which(x==1)[1])
    day.exposed[is.na(day.exposed)] <- 0
    
    day.expose.start <- apply(exposed.status,1, function(x) which(x==1)[1])
    day.expose.start[is.na(day.expose.start)] <- 0
    
    day.infect.start <- (day.expose.start + round(expose.dist))*(i>day.expose.start & (day.expose.start !=0 ))
    day.infect.end <- day.infect.start + round(infect.dist)*(i>day.expose.start & (day.expose.start !=0 ))
    
    n.infect.prev[,i] <- sum(infect.status[,(i-1)]) #how many people in HH were infectious at previous time?
    
    infect.status[,i]  <- (i >= day.infect.start) * (i <=day.infect.end ) * (day.infect.start>0)  #You are infectious for days in specified range
    
    exposed.status[,i] <- (1-rbinom(nrow(df1), 1, prob.uninf.day.comm*prob.uninfect.day^n.infect.prev[,i] )) ^ (1- exposed.status[,(i-1)]) #exponent ensure once you are exposed, you stay in that category
    
  } 
  
  ##Assume that get PCR on day 1 of being infectious #(2 days asymptomatic transmission)
  if(sum(day.expose.start)>0){
    earliest.expose.hh <- min(day.expose.start[day.expose.start>0])
  } else{
    earliest.expose.hh <- 1
  }
  latest.expose.hh <- min((earliest.expose.hh+20),ncol(exposed.status))
  
  start <- as.Date("2020-08-29")
  df1$PCR_DATE <- as.Date(day.infect.start,origin=start)+runif(1,min=1,max=100)
  df1$t.index <- as.numeric(df1$PCR_DATE - min(df1$PCR_DATE)) 
  df1$infected <- apply(exposed.status[,earliest.expose.hh:latest.expose.hh, drop=F],1,max) #Only count infections if they were exposed within 21 days of first exposure in HH
  return(df1)
}






