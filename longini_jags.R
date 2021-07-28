longini_jags  <- "
model{
  
  for(k in 1:max.hh.size){
    ds[1:(k+1),k] ~ dmulti(m[1:(k+1), k], n.hh[k])
  }
  
  #Could define B=1-P1; Q=1-P2; then have P1 and P2 be a function of covariates e.g., proportion vaccinated)
  logit_B ~ dnorm(0, 1e-4)
  logit_Q ~ dnorm(0, 1e-4)

  B <- exp(logit_B)/(1+ exp(logit_B))
  Q <- exp(logit_Q)/(1+ exp(logit_Q))

  
  m[1,1] = B   # Probability of 0 of 1 HH member infected
  m[2,1] = 1-B # Probability of 1 of 1 HH member infected

  for (k in 2:max.hh.size){
    m[1,k] = B^k # Probability everyone in HH escapes infection from the community 
    
    for (j in 1:(k-1)){
        m[j+1,k] = choose_kj_mat[k,j]*m[j+1,j]*(B^(k-j))*Q^(j*(k-j)) # Probabilty j out of k HH members infected
    }
  
     m[k+1,k] = 1-sum(m[1:k,k]) #Probability that everyone in HH infected
  }
}
"
