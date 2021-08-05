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
      countsUninf <-  d1$exposed.date > x 
      countsInf <-  (d1$exposed.date == x)*d1$infected  #how many people exposed this date? 
      
      grpN <- aggregate( cbind.data.frame(countsUninf,countsInf), 
                                     by = list(date=rep(x, length(vax)),vax1 = vax, agec=d1$agegrp), 
                                     FUN = sum)
      return(grpN)
    }
    ) 
      
    n.date.uninf <- bind_rows(n.date.uninf)
    
    #This plot shows the number of people in each category by date
    plot(n.date.uninf$date, n.date.uninf$countsUninf)
    plot(n.date.uninf$date, n.date.uninf$countsInf)
    
    ###NOTE THIS NEEDS TO BE MODIFIED TO COUNT NUMBER OF PEOPLE unINFECTED
    #AT EACH TIME POINT AND NUMBER OF PEOPLE INFEXTED AT EACH POINT

    
####################################################################################################################
#Step 2, Count exposure days in the household based on the number of infected people in each group
####################################################################################################################
    
    #NOTE: need to filter out HH with no infections
    
    #For each HH that has an infection, need to define for each day how many contacts are infected
    
    #What is the first date of infection in each household
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
    
    
    
 
    #could vectorize this whole thing instead of loop
   # for(d in 0:40){
    for(d in 3){
    
      b1$cal_date <- b1$hh.first.date + d #what is the current date in the HH (defined relative to first infection in HH)
      
      b1$vax <- ((b1$vax2dose_date + 10) < b1$cal_date )#was the person vaccinated 10+ days before the current date?
      
      b1$not_infected <- b1$exposed.date > b1$cal_date  # has the person been infected by cal_date?
      
      b1$infectious <- (b1$cal_date >= b1$infect.date) & (b1$cal_date < b1$end.infectious.date) #is the individual infectious at this time?
    
      c1 <- b1[,c("HH_CERTAIN",'agegrp', 'vax', 'not_infected', 'infectious')]
    
      #(Next step--copy c1, change names, then do man-many merge by HH with c1)
    }
    
    
