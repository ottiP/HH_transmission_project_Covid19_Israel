library(data.table)
library(dplyr)
d1 <- readRDS('./Data/simulated_data.rds')



start.date <- as.Date('2020-06-01')
end.date <- as.Date('2021-04-01')
all.date <- seq.Date(from=start.date, to=end.date, by='day')

#When exposed? (this is date of censoring)
d1$exposed.date <- d1$PCR_DATE - rgamma(nrow(d1), 7,2)
d1$exposed.date[d1$infected==0] <- end.date #if uninfected; censor at the end of follow up time

#When infectious?
d1$infect.date <- d1$exposed.date + rgamma(nrow(d1), 3,2)
d1$infect.date[d1$infected==0] <- end.date #if uninfected; censor at the end of follow up time

#When is the person no longer infectious
d1$end.infectious.date <- d1$infect.date + rgamma(nrow(d1), 7,1)

d1$agegrp <- 1
d1$agegrp[d1$AGE>18 & d1$age<60] <- 2
d1$agegrp[d1$AGE>=60 & d1$age<150] <- 3

####################################################################################################################
#Step 1, count all people who have not yet been infected, regardless of household infection status, 
#by vaccine status and age group at that time. This is used for estimating exogenous piece of the likelihood
####################################################################################################################

    n.date.uninf <- lapply(all.date, function(x){
      vax <- ((d1$vax2dose_date + 10) < x )#was the person vaccinated 10+ days before the current date?
      counts <- data.table( d1$exposed.date > x ) 
      grpN <- aggregate(x = counts, 
                                     by = list(date=rep(x, length(vax)),vax1 = vax, agec=d1$agegrp), 
                                     FUN = length)
      return(grpN)
    }
    ) 
      
    n.date.uninf <- bind_rows(n.date.uninf)
    
    #This plot shows the number of people in each category by date
    plot(n.date.uninf$date, n.date.uninf$V1)

    
####################################################################################################################
#Step 2, Count exposure days in the household based on the number of infected people in each group
####################################################################################################################
    

    #What is the first date of infection in each household, how many Infectins per HH, how many people per HH?
    d1 <- d1 %>%
      group_by(HH_CERTAIN) %>%
      mutate(
        hh.first.date = min(exposed.date, na.rm = T),
        n.infections.hh = sum(infected, na.rm = T),
        n.people.hh = length(infected),
        
      ) %>%
      arrange(HH_CERTAIN)
    
    #FIlter out HH with at least 1 infections and 2+ people for this calculation since focus is on HH transmission specifically
    b1 <- d1[d1$n.infections.hh>0 & d1$n.people.hh>=2 ,]
    
    
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
      e1$hh_risk_time <- e1$hh_risk_time1 * (e1$hh_risk_time1<e1$hh_risk_time2) + e1$hh_risk_time2 * (e1$hh_risk_time2<e1$hh_risk_time1)
      e1$hh_risk_time[e1$hh_risk_time <0] <- 0
      
      e1 <- e1[,c('agegrp','contactagegrp', 'vax','contactvax','hh_risk_time')]
      
      #then aggregate the hh_risk time by group
      
      f1 <- aggregate( e1[,'hh_risk_time'], by=list('agegrp'=e1$agegrp,'contactagegrp'=e1$contactagegrp, 'vax'=e1$vax,'contactvax'=e1$contactvax), FUN=sum)
      
      saveRDS(f1,'./Data/discrete_data.rds')
      