library(data.table)
library(dplyr)
d1 <- readRDS('./Data/simulated_data.rds')



start.date <- as.Date('2020-06-01')
end.date <- as.Date('2021-04-01')
all.date <- seq.Date(from=start.date, to=end.date, by='day')

d1$infect.date <- d1$PCR_DATE - rgamma(nrow(d1), 7,2)
d1$infect.date[d1$infected==0] <- end.date #if uninfected; censor at the end of follow up time

d1$agegrp <- 1
d1$agegrp[d1$AGE>18 & d1$age<60] <- 2
d1$agegrp[d1$AGE>=60 & d1$age<150] <- 3

#Step 1a, count all people who have not yet been infected, regardless of household infection status, 
#by vaccine status and age group at that time

n.date.uninf <- lapply(all.date, function(x){
  vax <- ((d1$vax2dose_date + 10) < x )#was the person vaccinated 10+ days before the current date?
  counts <- data.table( d1$infect.date > x ) 
  grpN <- aggregate(x = counts, 
                                 by = list(date=rep(x, length(vax)),vax1 = vax, agec=d1$agegrp), 
                                 FUN = length)
  return(grpN)
}
) 
  
n.date.uninf <- bind_rows(n.date.uninf)

#This plot shows the number of people in each category by date
plot(n.date.uninf$date, n.date.uninf$V1)
