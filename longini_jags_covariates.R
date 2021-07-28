longini_jags_covar  <- "
model{
  
  for(v in 1:2){ #vaccine grup
  
  for(k in 1:max.hh.size){
    ds[1:(k+1),k,v] ~ dmulti(m[1:(k+1), k,v], n.hh[k,v])
  }
  
  #Could define B=1-P1; Q=1-P2; then have P1 and P2 be a function of covariates e.g., proportion vaccinated)
  B[v] <- 1 - P1[v]
  Q[v] <- 1 - P2[v]
  
  logit_P1[v] <- alpha1 + (v==2)*beta1
  logit_P2[v] <- alpha2 + (v==2)*delta1
  
  P1[v] <- exp(logit_P1[v])/(1+ exp(logit_P1[v]))
  P2[v] <- exp(logit_P2[v])/(1+ exp(logit_P2[v]))
  
  m[1,1,v] = B[v]   # Probability of 0 of 1 HH member infected 
  m[2,1,v] = 1-B[v] # Probability of 1 of 1 HH member infected 
    
  for (k in 2:max.hh.size){
    m[1,k,v] = B[v]^k # Probability everyone in HH escapes infection from the community 
    
    for (j in 1:(k-1)){
        m[j+1,k,v] = choose_kj_mat[k,j]*m[j+1,j, v]*(B[v]^(k-j))*Q[v]^(j*(k-j)) # Probabilty j out of k HH members infected
    }
  
     m[k+1,k,v] = 1-sum(m[1:k,k,v]) #Probability that everyone in HH infected
  }
  }
  
  #Note prior here is a bit tighter than usual to prevent extreme values
  alpha1 ~ dnorm(0, 1/4)
  alpha2 ~ dnorm(0, 1/4)
  
  #Prior on alpha and beta
  beta1 ~ dnorm(0, 1/4) #(inv.var-0.04 == sd=5; middle%ci roughly ranges from -3.3,-3.3; which is on scale of ratio is [0.04, 27]
  delta1 ~ dnorm(0, 1/4)

}
"
