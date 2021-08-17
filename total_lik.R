total_lik <- function(params,Y,X,cases,df_noinf,cases1,s0,s1,s2,s0_sum,s1_sum,s2_sum,vax_cats,n_vax_cats){
  ll_inf <- chain_bin_lik(params,Y,X)
  ll_noinf <- lik_noinf(params,cases,df_noinf)
  ll_singleinf<-lik_singleHH_inf(cases1,s0,s1,s2,s0_sum,s1_sum,s2_sum,vax_cats,n_vax_cats)
  ll_total <- ll_inf+ll_noinf+ll_singleinf
  return(ll_total)
}