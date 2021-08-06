#################################################################################################################################
#Likelihood Combinations and Calculations: define function
#################################################################################################################################
lik_noinf <- function(params,cases,df_noinf){
  
  library(lubridate)
  df_noinf$SECOND_VACCINE_DATE[is.na(df_noinf$vax2dose_date) == 1]<-enddate + 1  #Unvaccinated people (day after study ends)
  df_noint$SECOND_VACCINE_DATE<-floor_date(df_noint$SECOND_VACCINE_DATE)
  vax_cats<-unique(df_noinf$vax2dose_date)
  n_vax_cats<-length(vax_cats)
  
  start<-floor_date(start)
  enddate<-floor_date(enddate)
  
  age1<-as.numeric(df_noinf$AGE >= 10 & df_noinf$AGE < 60)
  age2<-as.numeric(df_noinf$AGE >= 60)
  
  age0_log_probs<-rep(0.00,
                      times = n_vax_cats)
  age1_log_probs<-rep(0.00,
                      times = n_vax_cats)
  age2_log_probs<-rep(0.00,
                      times = n_vax_cats)
  for(j in 1:n_vax_cats){
    print(j)
    #age0
    x<-cbind(0,
             1,  
             c(rep(0, times = as.numeric(vax_cats[j] - start)), rep(1, times = as.numeric(enddate - vax_cats[j] + 1))),
             0,
             0,
             0,
             cases)
    p0<-1.00/(1.00 + exp(-x%*%params))
    age0_log_probs[j]<-sum((df_noinf$vax2dose_date == vax_cats[j]) & (age1 == 0) & (age2 == 0))*sum(log(1.00 - p0))
    
    #age1
    x<-cbind(0,
             1,
             c(rep(0, times = as.numeric(vax_cats[j] - start)), rep(1, times = as.numeric(enddate - vax_cats[j] + 1))),
             0,
             1,
             0,
             cases)
    p0<-1.00/(1.00 + exp(-x%*%params))
    age1_log_probs[j]<-sum((df_noinf$vax2dose_date == vax_cats[j]) & (age1 == 1))*sum(log(1.00 - p0))
    
    #age2
    x<-cbind(0,
             1,
             c(rep(0, times = as.numeric(vax_cats[j] - start)), rep(1, times = as.numeric(enddate - vax_cats[j] + 1))),
             0,
             0,
             1,
             cases)
    p0<-1.00/(1.00 + exp(-x%*%params))
    age2_log_probs[j]<-sum((df_noinf$vax2dose_date == vax_cats[j]) & (age2 == 1))*sum(log(1.00 - p0))
    
  }
  ll_contribution<-sum(age0_log_probs +
                         age1_log_probs +
                         age2_log_probs)
  
  return(-ll_contribution)
  
}







