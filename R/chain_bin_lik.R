
#Note here, we have all time points represented in the df, so the likelihood is very simple--no need to exponentiate stuff  
chain_bin_lik_agg <- function(par,first_grp,X_hh,X_comm,group_id,new.inf,hh_par_index,comm_par_index, N_susceptible, N_inf_contacts){
  
  #### Define logit_p = X*params; need  to add as.matrix
  log_p_hh <- as.vector(as.matrix(X_hh) %*% par[hh_par_index]) ## Added as.matrix
  
  p_hh=  exp(log_p_hh) *(N_inf_contacts>0) + 0.5*(N_inf_contacts==0) #set to 0.5 if N.inf.contacts=0 to avoid log(0)...this doesn't influence results at all
  
  #Community pieces
  log_p_comm <- as.vector(as.matrix(X_comm) %*% par[comm_par_index]) ## Added as.matrix
  
  p_comm <- exp(log_p_comm)
  
  #If a person has multiple contacts in different risk groups, this will be spread over multiple rows..p_comm should only count first
  log.qi.row <-    first_grp*log(1-p_comm) + N_inf_contacts*log(1-p_hh)
  
  df.row <- cbind.data.frame(N_susceptible, log.qi.row, new.inf,group_id )
  
  suppressMessages( df.person <- df.row %>%
                      group_by(group_id, new.inf, N_susceptible) %>%
                      summarize(log.qi.person = sum(log.qi.row))
  ) 
  #sum(df.person$N_susceptible[df.person$new.inf==1])
  #LogSumExp trick.  add log()
  log.pi.person <-  log(1 - exp(df.person$log.qi.person))
  
  log.ll.piece= df.person$N_susceptible*(df.person$new.inf*log.pi.person + (1-df.person$new.inf)*df.person$log.qi.person) 
  
  if(is.nan(sum(log.ll.piece ))){
    ll= -1e8
  }else{
    ll <- sum(log.ll.piece )
  }
  # }
  
  return(-ll)
}
