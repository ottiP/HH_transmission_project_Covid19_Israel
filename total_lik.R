total_lik <- function(params,Y,X,df_noinf,cases){
  ll_inf <- chain_bin_lik(params,Y,X)
  ll_noinf <- lik_noinf(params,cases,df_noinf)
  ll_total <- ll_inf+ll_noinf
  return(ll_total)
}