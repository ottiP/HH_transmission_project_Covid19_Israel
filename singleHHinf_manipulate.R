######################################################
#Fake Data
######################################################
singleHHinf_manipulate <- function(df){
  complete_dates<-rep(startdate,
                      times = (enddate - startdate + 1))
  for(j in 1:(enddate - startdate + 1)){
    complete_dates[j]<-as.Date(startdate + j - 1)
  }
  
  #################################################################################################################################
  #Likelihood Combinations and Calculations
  #################################################################################################################################
  df$SECOND_VACCINE_DATE <- ifelse((df$SECOND_VACCINE_DATE+10)<=as.character(enddate),as.character(df$SECOND_VACCINE_DATE+10),as.character(enddate+1))
  df$SECOND_VACCINE_DATE[is.na(df$SECOND_VACCINE_DATE) == 1]<-as.character(enddate + 1)  #Unvaccinated people (day after study ends)
  df$SECOND_VACCINE_DATE<-floor_date(as.Date(df$SECOND_VACCINE_DATE))
  vax_cats<-unique(df$SECOND_VACCINE_DATE)
  n_vax_cats<-length(vax_cats)
  
  start<-floor_date(startdate)
  enddate<- floor_date(enddate)
  complete_dates<-floor_date(complete_dates)
  
  age1<-as.numeric(df$AGE >= 12 & df$AGE < 60)  #CHECK AGE GROUPS
  age2<-as.numeric(df$AGE >= 60)                            #CHECK AGE GROUPS
  
  df$exposed.date <- floor_date(df$PCR_DATE - rgamma(nrow(df),shape=4,scale=1)-rgamma(nrow(df), shape=2.34,scale=1/(2.59)), unit='day') #note: need to use floor date otherwise have fractional days, which causes problems later
  ###############################################################################################################################################################
  #Getting Counts:
  #This should be calculated outside of the likelihood function (data preparation stage) and "s0, s1, s2, s0_sum, s1_sum, s2_sum" should be input to the function
  ###############################################################################################################################################################
  s0<-matrix(0.00,
             nrow = n_vax_cats,
             ncol = (enddate - start + 1))
  s1<-matrix(0.00,
             nrow = n_vax_cats,
             ncol = (enddate - start + 1))
  s2<-matrix(0.00,
             nrow = n_vax_cats,
             ncol = (enddate - start + 1))
  s0_sum<-matrix(0.00,
                 nrow = n_vax_cats,
                 ncol = (enddate - start + 1))
  s1_sum<-matrix(0.00,
                 nrow = n_vax_cats,
                 ncol = (enddate - start + 1))
  s2_sum<-matrix(0.00,
                 nrow = n_vax_cats,
                 ncol = (enddate - start + 1))
  for(j in 1:n_vax_cats){
    
    for(k in 1:(enddate - start + 1)){
      
      s0[j,k]<-sum((df$SECOND_VACCINE_DATE == vax_cats[j]) & (age1 == 0) & (age2 == 0) & (round(df$exposed.date) == round(complete_dates[k])))  #SHOULD BE INFECTION DATE (AFTER DELAY DISTN ADJUSTMENT)
      s1[j,k]<-sum((df$SECOND_VACCINE_DATE == vax_cats[j]) & (age1 == 1) & (round(df$exposed.date) == round(complete_dates[k])))                #SHOULD BE INFECTION DATE (AFTER DELAY DISTN ADJUSTMENT)
      s2[j,k]<-sum((df$SECOND_VACCINE_DATE == vax_cats[j]) & (age2 == 1) & (round(df$exposed.date) == round(complete_dates[k])))                #SHOULD BE INFECTION DATE (AFTER DELAY DISTN ADJUSTMENT)
      
    }
  }
  
  for(j in 1:n_vax_cats){
    for(k in 1:(enddate- start)){
      
      s0_sum[j,k]<-sum(s0[j, (k + 1):(enddate - start + 1)])
      s1_sum[j,k]<-sum(s1[j, (k + 1):(enddate - start + 1)])
      s2_sum[j,k]<-sum(s2[j, (k + 1):(enddate - start + 1)])
      
    }
  }
  
  out.list <- list('s0'=s0,'s1'=s1,'s2'=s2,'s0_sum'=s0_sum,'s1_sum'=s1_sum,'s2_sum'=s2_sum,'vax_cats'=vax_cats,'n_vax_cats'=n_vax_cats)
  return(out.list)
}


