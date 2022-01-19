#' Title
#'
#' @param N 
#' @param hh.size 
#' @param p.comm.base 
#' @param p.hh.base 
#' @param prop.vax 
#' @param delay 
#' @param VE_susceptibility 
#' @param VE_txn 
#' @param duration.latent 
#' @param duration_infectious 
#' @param p.detect 
#'
#' @return
#' @export
#'
#' @examples
sim.hh.func.misclassified <- function(N,
                                      #n.hh = 1000,
                                      hh.size = 3,
                                      p.comm.base = 0.0002,
                                      p.hh.base = 0.02,
                                      prop.vax = 0.5,
                                      delay=T,
                                      
                                      VE_susceptibility = 90, #vaccine effectiveness against infection
                                      VE_txn = 50, #VE infectiousness
                                      
                                      duration.latent = 4,
                                      duration_infectious = 5,
                                      #What is probability that an infected person is detected on each day 
                                      p.detect = 0.99 #what is overall probability that an infection is detected?
){
  
  p.detect.day <- 0.75#1 - exp(log(1-p.detect)/duration_infectious) 
  
  time.steps = 300
  
  latent <- array(0, dim=c(4,time.steps,hh.size))
  
  infectious <- array(0, dim=c(4,time.steps,hh.size))
  
  immune <- matrix(NA, nrow=time.steps, ncol=hh.size)
  #vax <- matrix(NA, nrow=time.steps, ncol=hh.size)
  
  
  vax.status <- rbinom(hh.size,1,prop.vax) #is individual vaccinated?--note this is not time varying
  
  
  latent[1,1,] <- rbinom(hh.size, 1, p.comm.base)#1st of the 4 latent classes
  
  infectious[1,1,] <- 0 #first of the 4 infectious classes
  
  immune[1,] <- 0
  
  #vax[1,] <- 0
  #Contribution of community infection
  log.p.comm <- log(p.comm.base) + vax.status*log((1-VE_susceptibility/100))
  p.comm <- exp(log.p.comm)
  
  #Contribution of HH infection from unvaccinated contact, per contact
  log.p.hh.unvax <- log(p.hh.base) + vax.status*log((1-VE_susceptibility/100))
  p.hh.unvax <- exp(log.p.hh.unvax) 
  
  #Contribution of HH infection from vaccinated contact, per contact
  log.p.hh.vax <- log(p.hh.base) + vax.status*log((1-VE_susceptibility/100)) +
    log(1-VE_txn/100)
  p.hh.vax <- exp(log.p.hh.vax)
  
  for(i in 2:time.steps){
    for(j in 1:hh.size){
      
      #need to exponentiate these by the number of infectious HH members in each category
      N.infectious.contact <- sum(infectious[,(i-1),-j])
      N.vaxed.infected.contact <- sum(infectious[,(i-1),-j] %*% diag(vax.status[-j]))
      N.unvax.infectious.contact <- N.infectious.contact - N.vaxed.infected.contact
      
      log.qi <- log(1 - p.comm[j]) + N.vaxed.infected.contact*log(1 - p.hh.vax[j]) + N.unvax.infectious.contact*log(1-p.hh.unvax[j])
      
      qi <- exp(log.qi)
      
      pi <- 1- qi #probability of being infected at the time point
      
      susceptible <- (1-sum(infectious[,(i-1),j])) * (1-immune[(i-1),j]) * (1-sum(latent[,(i-1),j])) #is person susceptible at previous time step
      
      new.inf <-  rbinom(1,1,pi)*susceptible #is a person who was previously susceptibly now latent?
      
      #If a person has been latent for <4 time points, stay latent in current time point
      # stay.latent <-   (sum(latent[1:(i-1),j])<4) 
      # 
      # #If a person has been infectious for <5 time points, stay infectious
      # stay.infectious <- (sum(infectious[1:(i-1),j])<5) 
      # 
      
      #Can use same value regardless of whether person was in latent 1,2 3 etc
      stay.latent <-   sum(latent[,(i-1),j])*(rbinom(1,1, (1-1/duration.latent)^1/4))  #Transition to infectious if latent previously
      
      stay.infectious <- sum(infectious[,(i-1),j])*(rbinom(1,1, (1-1/duration_infectious)^1/4))
      
      
      #transition between latent states
      latent[1,i,j] <- (new.inf  + #if susceptible at previous time step, do they become infected (latent)
                          latent[1,(i-1),j]*stay.latent ) #If latent previously, do they stay latent?
      
      latent[2,i,j] <- latent[1,(i-1),j]*(1-stay.latent ) + #If latent previously, do they stay latent?
        latent[2,(i-1),j]*stay.latent 
      
      latent[3,i,j] <- latent[2,(i-1),j]*(1-stay.latent ) + #If latent previously, do they stay latent?
        latent[3,(i-1),j]*stay.latent 
      
      latent[4,i,j] <- latent[3,(i-1),j]*(1-stay.latent ) + #If latent previously, do they stay latent?
        latent[4,(i-1),j]*stay.latent 
      
      infectious[1,i,j]   <-     latent[4,(i-1),j]*(1-stay.latent) +  #if latent at previous time step, do they become infectious?
        infectious[1,(i-1),j]*(stay.infectious)     #if infectious at previous step, do they become immune
      
      infectious[2,i,j]   <-    infectious[1,(i-1),j]*(1-stay.infectious) +  #if latent at previous time step, do they become infectious?
        infectious[2,(i-1),j]*(stay.infectious)     #if infectious at previous step, do they become immune
      
      infectious[3,i,j]   <-    infectious[2,(i-1),j]*(1-stay.infectious) +  #if latent at previous time step, do they become infectious?
        infectious[3,(i-1),j]*(stay.infectious)
      
      infectious[4,i,j]   <-    infectious[3,(i-1),j]*(1-stay.infectious) +  #if latent at previous time step, do they become infectious?
        infectious[4,(i-1),j]*(stay.infectious)
      
      
      immune[i,j]   <-     infectious[4,(i-1),j] *(1-stay.infectious) + #if previously infectious, do they become immune?
        immune[i-1,j] #if previously immune, stay immune
      #vax[i,j]<- ifelse(vax.status[j]==1,1,0)
    }
  }
  
  #Observation probability
  detect.inf <- infectious * rbinom(length(infectious),1, p.detect.day ) #detection of case ?
  # if(delay==T){
  #   detect.inf <- infectious
  # }else{
  #   detect.inf <- latent
  # 
  # }
  
  if(hh.size>1){
    first.infection.date <- apply(detect.inf[1,,], 2, function(x){
      if(sum(x)>0){
        z <- which(x==1)[1] 
      }else{
        z=-9999
      }
      return(z)
    })
  }else{
    first.infection.date <- apply(detect.inf[1,,,drop=F], c(1,3), function(x){
      if(sum(x)>0){
        z <- which(x==1)[1] 
      }else{
        z=-9999
      }
      return(z)
    })
  }
  
  true_duration_infectious <- apply(infectious, c(3), function(x){
    if(sum(x)>0){
      z = sum(x)  
    }else{
      z=-9999
    }
    return(z)
  })
  
  out.df <- cbind.data.frame(vax.status,first.infection.date,true_duration_infectious, 'HH'=N)
  out.df$infected <- 1*( out.df$first.infection.date != -9999)
  
  return(out.df)
  
  #out.ls <- list('latent'=latent,'infectious'=infectious,'immune'=immune ,'first.infection.date'=first.infection.date,'vax.status'=vax.status)
  
}
