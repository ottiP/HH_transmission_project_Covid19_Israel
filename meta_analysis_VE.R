meta_analysis_func  <- function(VE,SE2,N_iter){
  
  model_string<- "
   model{
    for(i in 1:N_iter){
      VE_hat[i] ~ dnorm(VE_iter[i],inv_var[i])
      VE_iter[i] ~ dnorm(VE,tau2)
    }
    
    VE ~ dnorm(0,0.0001)
    tau2 = 1/(sigma*sigma)
    sigma ~ dunif(0,1000)
    
   }"
  
  #### #### #### #### #### #### #### #### #### #### #### #### 
  #### Model fitting
  #### #### #### #### #### #### #### #### #### #### #### #### 
  
  inits1 =  list(".RNG.seed"=c(123),".RNG.name"='base::Wichmann-Hill')
  inits2 =  list(".RNG.seed"=c(456),".RNG.name"='base::Wichmann-Hill')
  inits3 =  list(".RNG.seed"=c(789),".RNG.name"='base::Wichmann-Hill')
  
  model_spec<- textConnection(model_string)
  model_jags <- jags.model(model_spec,
                           inits=list(inits1,inits2,inits3),
                           data=list('VE_hat'=VE,
                                     'inv_var'=1/(SE2),
                                     'N_iter'=N_iter),
                           n.adapt=10000,
                           n.chains=3)
  params<-c('VE',
            'tau2')
  
  
  #### #### #### #### #### #### #### #### #### #### #### #### 
  #### Posterior Sampling
  #### #### #### #### #### #### #### #### #### #### #### #### 
  posterior_samples  <- coda.samples(model_jags, 
                                     params,
                                     n.iter=100000)
  posterior_samples.all <- do.call(rbind,posterior_samples)
  
  #### #### #### #### #### #### #### #### #### #### #### #### 
  #### Posterior Inference
  #### #### #### #### #### #### #### #### #### #### #### #### 
  post_means <- apply(posterior_samples.all,2,median)
  par.names <- names(post_means)
  ci<-t(hdi(posterior_samples.all,credMass=0.95))
  ci<-matrix(sprintf("%.1f",round(ci,1)),ncol=2)
  row.names(ci) <- par.names
  post_means<-sprintf("%.1f",round(post_means,1))
  names(post_means) <- par.names
  VE<-c(post_means[1],ci[1,])
  tau2<-c(post_means[2],ci[2,])
  res.list <- list("VE"=VE,"tau2"=tau2)
  
  return(res.list)
}





weights_VE <- function(ds){
  
  # Standardized VE against susceptibility
  f_u <- (sum(ds$m2[ds$Y==1]))/(sum(ds$Y))
  f_v <- (sum(ds$z2[ds$Y==1]))/(sum(ds$Y))
  f_e <- 1-f_v-f_u
  
  # Standardized VE against infectiousness
  ds_vax <- ds[ds$vax2a==1,]
  g_v = (sum(ds_vax$Y[(ds_vax$m2==1)|(ds_vax$z2==1)]))/(sum(ds$Y[(ds$m2==1)|(ds$z2==1)]))
  g_u <- 1-g_v
  out.list <- list("f_u"=f_u,"f_v"=f_v,"f_e"=f_e,"g_u"=g_u,"g_v"=g_v)
  return(out.list)
}













