generate_delay_dist_exp<-function(df){
  n.sim <- nrow(df)
  #Could modify one of these to be time from infection to test
  expose.dist= rexp(n.sim, 1/4) #duration latent: Davies et al., The Lancet, 2020
  infect.dist= rexp(n.sim, 1/5) #duration of subclinical infectiousness -> temptative, this should be duration of clinical infectiousness, Davies, Lancet, 2020
  delay.dist= rexp(n.sim, 1/0.9) #duration infectiousness to test  #from Li.et al., Science  May 2020, estimates after 23rd Jan, parameters for Gamma should be a and b
  
  df$date.exposed <- floor_date(as.Date(df$PCR_DATE - delay.dist - expose.dist),'day')
  df$date.exposed[df$IS_POSITIVE_CD!=2] <- as.Date('2099-01-01')
  
  df$date.infected.end <-  floor_date(as.Date(df$PCR_DATE + infect.dist - delay.dist),'day')
  df$date.infected.end[df$IS_POSITIVE_CD!=2] <- as.Date('2099-01-01')
  
  df$date.onset <-  floor_date(as.Date(df$PCR_DATE - delay.dist),'day')
  df$date.onset[df$IS_POSITIVE_CD!=2] <- as.Date('2099-01-01')
  
  df$id <- 1:nrow(df)
  df$vax <- df$vax.status
  
  df <- df[, c("date.exposed" ,"IS_POSITIVE_CD", "date.infected.end", "date.onset" ,"id","HH",'vax')] 
  
  df <- df %>%
    group_by(HH) %>%
    mutate(N_people_HH= n_distinct(id), N_infect_hh= sum(IS_POSITIVE_CD==2))
  
  return(df)
}
