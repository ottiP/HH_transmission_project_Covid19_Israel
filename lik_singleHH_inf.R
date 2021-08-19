####################################################################
#Reading in Data
####################################################################

###############################################################################################################################################################
# Likelihood function for single HH with infections 
lik_singleHH_inf <- function(cases1,s0,s1,s2,s0_sum,s1_sum,s2_sum,vax_cats,n_vax_cats){
  
  start<- floor_date(startdate)
  enddate<- floor_date(enddate)
  age0_log_probs<-rep(0.00,
                      times = n_vax_cats)
  age1_log_probs<-rep(0.00,
                      times = n_vax_cats)
  age2_log_probs<-rep(0.00,
                      times = n_vax_cats)
  for(j in 1:n_vax_cats){
    #print(j)
    #age0
    x<-cbind(0,
             1,  
             c(rep(0, times = as.numeric(vax_cats[j] - start)), rep(1, times = as.numeric(enddate - vax_cats[j] + 1))),
             0,
             0,
             0,
             cases1)
    p0<-1.00/(1.00 + exp(-x%*%params))
    age0_log_probs[j]<-sum(s0[j,]*log(p0)) +
      sum(s0_sum[j,]*log(1.00 - p0))
    
    #age1
    x<-cbind(0,
             1,
             c(rep(0, times = as.numeric(vax_cats[j] - start)), rep(1, times = as.numeric(enddate - vax_cats[j] + 1))),
             0,
             1,
             0,
             cases1)
    p0<-1.00/(1.00 + exp(-x%*%params))
    age1_log_probs[j]<-sum(s1[j,]*log(p0)) +
      sum(s1_sum[j,]*log(1.00 - p0))
    
    #age2
    x<-cbind(0,
             1,
             c(rep(0, times = as.numeric(vax_cats[j] - start)), rep(1, times = as.numeric(enddate - vax_cats[j] + 1))),
             0,
             0,
             1,
             cases1)
    p0<-1.00/(1.00 + exp(-x%*%params))
    age2_log_probs[j]<-sum(s2[j,]*log(p0)) +
      sum(s2_sum[j,]*log(1.00 - p0))
    
  }
  ll_contribution<-sum(age0_log_probs +
                         age1_log_probs +
                         age2_log_probs)
  
  
  return(-ll_contribution)
  
}


















