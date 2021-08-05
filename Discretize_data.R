library(data.table)
library(dplyr)
library(lubridate)
a1 <- readRDS('./Data/simulated_data.rds')
set.seed(123)


start.date <- as.Date('2020-06-01')
end.date <- as.Date('2021-04-01')
all.date <- seq.Date(from=start.date, to=end.date, by='day')

#When exposed? (this is date of censoring)
a1$exposed.date <- floor_date(a1$PCR_DATE - rgamma(nrow(a1), 7,2), unit='day') #note: need to use floor date otherwise have fractional days, which causes problems later
a1$exposed.date[a1$infected==0] <- end.date #if uninfected; censor at the end of follow up time

#When infectious?
a1$infect.date <- floor_date(a1$exposed.date + rgamma(nrow(a1), 3,2), unit='day')
a1$infect.date[a1$infected==0] <- end.date #if uninfected; censor at the end of follow up time

#When is the person no longer infectious
a1$end.infectious.date <- floor_date(a1$infect.date + rgamma(nrow(a1), 7,1), unit='day')

a1$agegrp <- 1
a1$agegrp[a1$AGE>18 & a1$age<60] <- 2
a1$agegrp[a1$AGE>=60 & a1$age<150] <- 3

#What is the first date of infection in each household, how many Infectins per HH, how many people per HH?
a1 <- a1 %>%
  group_by(HH_CERTAIN) %>%
  mutate(
    hh.first.date = min(exposed.date, na.rm = T),
    n.infections.hh = sum(infected, na.rm = T),
    n.people.hh = length(infected),
    
  ) %>%
  arrange(HH_CERTAIN)

a1$hh_infected <- a1$n.infections.hh > 1 #Is this a HH that is ultimately infected?

a1$hh_single <- a1$n.people.hh ==1 #is this a single person HH


####################################################################################################################
#Step 1, count all people who have not yet been infected, regardless of household infection status, 
#by vaccine status and age group at that time. This is used for estimating exogenous piece of the likelihood
####################################################################################################################

#Not sure if this is right--do we need to separate this out for households who ultimately become infected and those that don't?

n.date.uninf <- lapply(all.date, function(x){
  vax <- ((a1$vax2dose_date + 10) < x )#was the person vaccinated 10+ days before the current date?
  countsUninf <-  a1$exposed.date > x 
  countsInf <-  (a1$exposed.date == x)*a1$infected  #how many people exposed this date? 
  
  grpN <- aggregate( cbind.data.frame(countsUninf,countsInf), 
                     by = list(date=rep(x, length(vax)),vax1 = vax, agec=a1$agegrp, person_infected=a1$infected, hh_single=a1$hh_single), 
                     FUN = sum)
  return(grpN)
}
) 

n.date.uninf <- bind_rows(n.date.uninf)

n.date.uninf$cases.day <- rpois(nrow(n.date.uninf), lambda=100)  #NOTE: IN REAL DATA MERGE IN INCIDENCE DATA HERE

#This plot shows the number of people in each category by date
plot(n.date.uninf$date, n.date.uninf$countsUninf)
plot(n.date.uninf$date, n.date.uninf$countsInf)


## Calculate probability for each stratum

X.exo <- model.matrix(  ~ as.factor(agec) + as.factor(vax1) + cases.day  , data=n.date.uninf)

params.exo <- rep(-0.001, ncol(X.exo))  #Just plug in random values...this is part of what is estimated by mle

n.date.uninf$log_r_c <- as.vector(X.exo %*% params.exo)

n.date.uninf$q_exo <-   1 - exp(n.date.uninf$log_r_c)

#####################################################################################################
#LL for uninfected PERSON (rather than uninfected hh)
#########################################################################################################

uninf.hh <- n.date.uninf[n.date.uninf$person_infected==F,]

ll_exo_piece0 <- sum(  uninf.hh$countsUninf * log( ( uninf.hh$q_exo ))) ###LL contribution of uninfected PEOPLE ####


###########################################################################################
#Then LL piece for PEOPLE with infection 
###########################################################################################                             
inf.hh <- n.date.uninf[n.date.uninf$person_infected==T ,]
inf.hh.spl <- split(inf.hh, paste(inf.hh$vax1, inf.hh$agec))
inf.hh.spl <- lapply(inf.hh.spl, function(x) x[order(x$date),] )
inf.hh.spl <- lapply(inf.hh.spl, function(x){
  x$cum.q.t1 <-  exp(cumsum( log( x$q_exo) )) / x$q_exo #product of probability 1 -> t-1  
  x$prob.piece <- ((1- x$q_exo) * x$cum.q.t1  ) ^ x$countsInf 
  x <- x[x$countsInf >0,] #Only keeps strata/days when there is actually a case reported (other days used for calculating cum probabilities)
  return(x)
} 
)

inf.hh2 <- bind_rows(inf.hh.spl)

ll_exo_piece1 <- sum( log(inf.hh2$prob.piece))  ###LL contribution of infected HH ####

ll_exo <- ll_exo_piece0 + ll_exo_piece1  #Exogenous LL contribution




####################################################################################################################
#Step 2, Count exposure days in the household based on the number of infected people in each group
####################################################################################################################
    

     #FIlter out HH with no infections infections and 1 people for this calculation since focus is on HH transmission specifically
    b1 <- a1[a1$n.infections.hh>0 & a1$n.people.hh>=2 ,]
    
    
     ##ASSUMPTION: fix the vaccine status for each individual to its value at date of first infection in the HH (reasonable given relatively short time period)
      b1$vax <- ((b1$vax2dose_date + 10) < b1$hh.first.date )#was the person vaccinated 10+ days before the first infection in the HH?
      
      b1$vax[is.na(b1$vax)] <- 0 #if vaxdate=NA, vax=NA
      
      #b1$not_infected <- b1$exposed.date > b1$cal_date  # has the person been infected by cal_date?
      
     # b1$infectious <- (b1$cal_date >= b1$infect.date) & (b1$cal_date < b1$end.infectious.date) #is the individual infectious at this time?
    
      c1 <- b1[,c('ID',"HH_CERTAIN",'agegrp', 'vax', 'exposed.date', 'infect.date', 'end.infectious.date', 'infected')]
    
      d1 <- c1  #copies c1
      
      names(d1) <- paste0('contact',names(c1))
      
      #Note there is probably a more clever way to code this so we don't have to do many to many on all possible combos and then filter
      e1 <- merge( c1, d1, by.x='HH_CERTAIN', by.y='contactHH_CERTAIN') #many-to-many merge by HH
      
      e1 <- e1[e1$ID != e1$contactID,] #get rid of self-self rows
      
      e1 <- e1[e1$contactinfected == 1, ] #a pair can only contribute person time if 'contact' is infected at some point
    
      #Now count number of follow up days for the pairs
      e1$hh_risk_time1 <- (e1$contactend.infectious.date - e1$contactinfect.date)  #total time contact is infectious
      e1$hh_risk_time2 <-           ( e1$exposed.date - e1$contactinfect.date  )    ##time between when contact is infection and when person 1 is censored
    
      #Take the minimum of hh_risk_time1, hh_risk_time2
      e1$hh_risk_time <- as.numeric(e1$hh_risk_time1 * (e1$hh_risk_time1<e1$hh_risk_time2) + e1$hh_risk_time2 * (e1$hh_risk_time2<e1$hh_risk_time1))
     
       e1$hh_risk_time[e1$hh_risk_time <0] <- 0
      
      e1 <- e1[,c('agegrp','contactagegrp', 'vax','contactvax','hh_risk_time','infected')]
      
  
      #then aggregate the hh_risk time by group
      
      f1 <- aggregate( e1[,'hh_risk_time'], by=list('agegrp'=e1$agegrp,'contactagegrp'=e1$contactagegrp, 'vax'=e1$vax,'contactvax'=e1$contactvax,'infected'=e1$infected), FUN=sum)
      
      saveRDS(f1,'./Data/discrete_data.rds')
      write.csv(f1,'./Data/hh_follow_up.csv')
      write.csv(n.date.uninf,'./Data/exogenous_follow_up.csv')
      