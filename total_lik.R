total_lik <- function(params,Y,X,df1,df_noinf,cases,s0,s1,s2,s0_sum,s1_sum,s2_sum,first_day,last_day){
  ll_inf <- chain_bin_lik(params,Y,X)+lik_noinf(params,df1,first_day,last_day)
  ll_noinf <- lik_noinf(params,df_noinf,start,enddate)
  ll_singleinf<-lik_singleHH_inf(s0,s1,s2,s0_sum,s1_sum,s2_sum)
  ll_total <- ll_inf+ll_noinf+ll_singleinf
  return(ll_total)
}
