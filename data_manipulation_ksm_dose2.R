delay.gen.ksm.singleHH.age <- function(input_df){

  df1b  <-  suppressWarnings(delay_dist_sim_updated_age(input_df))
  
  
  df1b <- df1b[,c('ID','first.date.hh','HH_CERTAIN','vax1dose_date','vax2dose_date',
                  'infected','AGE','day.exposed','day.infectious','day.infectious.end','max.time.hh')]
  
  ### Create a  copy of df1 and call it df2
  df2 <- df1b
  
  ### Rename all of the variables on df2 as HH_CERTAIN_b, Person_ID_b, Vax_dose1b
  colnames(df2) <- paste(colnames(df2), 'b', sep='_')
  
  ### Merge df1 and df2 by HH_CERTAIN; name this df3 
  df3 <- merge(df1b,df2,by.y='HH_CERTAIN_b' , by.x='HH_CERTAIN')

  df3.spl <- split(df3, paste0(df3$ID, df3$ID_b))
  df4.spl <- lapply(df3.spl, function(x){
    x$exxpand.t <- NA
    x$expand.t[x$infected==0] <- x$max.time.hh[x$infected==0] #censor uninfected person at max.time.hh
    x$expand.t[x$infected==1] <- x$day.exposed[x$infected==1] #censor infected person at (latent) day of exposure
    
    x$rowN <- 1:nrow(x)
    df.t <- x[rep(x$rowN, times=x$expand.t),] #expands the df to number of time points
    
    #Add an index for each person--note these are all pegged to 
    #original t.index but shifted
    df.t$t.index <- unlist(lapply(x$expand.t, function(x){ 
      z <- 1:x 
      return(z)
    }))
    
    # keep rows where the time index is within range of the infectious period for person b;
    ### CHECK: 
    df.t$keep <-  T*ifelse(df.t$ID == df.t$ID_b,1,0)  | (df.t$t.index >= df.t$day.infectious_b & df.t$t.index <= df.t$day.infectious.end_b )*ifelse(df.t$infected_b==1,1,0)
    dft.k <- df.t[df.t$keep==1,]
    dft.k <- dft.k[,c('ID','ID_b','HH_CERTAIN','vax1dose_date','vax2dose_date','vax1dose_date_b','vax2dose_date_b',
                      'infected','AGE','infected_b','first.date.hh', 't.index','day.exposed')]
    return(dft.k)
  })
  
  df4 <-bind_rows(df4.spl)
  df4$date.current <- df4$first.date.hh + df4$t.index -1
  
  #Create elements for design matrix
  df4$delta0 <- 0
  df4$delta0[(df4$ID == df4$ID_b)] <- 1 #exogenous/community risk 
  
  df4$alpha0 <- 0
  df4$alpha0[df4$ID != df4$ID_b ] <- 1
  
  df4$vax1a <- (df4$date.current > (df4$vax1dose_date+10))*1
  df4$vax1a[is.na(df4$vax1dose_date)] <- 0
  # 
  # #Is contact/person b vaccinated at time?
  df4$vax1b <- (df4$date.current > (df4$vax1dose_date_b  +10))*1*df4$alpha0
  df4$vax1b[is.na(df4$vax1dose_date_b)] <- 0
  
  #Same for 2nd dose--is person A and/or person B vaccinated?
  df4$vax2a <- (df4$date.current > (df4$vax2dose_date  +10))*1
  df4$vax2a[is.na(df4$vax2dose_date)] <- 0

  #Is contact/person b vaccinated at time?
  df4$vax2b <- (df4$date.current > (df4$vax2dose_date_b  +10))*1*df4$alpha0
  df4$vax2b[is.na(df4$vax2dose_date_b)] <- 0  
  
  # df4$AGE[(df4$AGE<=10)]<- 0
  df4$AGE1 <- 0
  df4$AGE1[(df4$AGE<60)&(df4$AGE>10)]<- 1
  df4$AGE2 <- 0
  df4$AGE2[(df4$AGE>=60)]<- 1
  
  
  df4$infect.at.timet <- df4$infected
  df4$infect.at.timet[df4$t.index < df4$day.exposed] <- 0
  df4 <- df4[!is.na(df4$HH_CERTAIN),]
  
  df4$time_risk <- as.numeric(pos_test_complete[df4$t.index])*df4$delta0
  df4$time_risk <- scale(df4$time_risk)
  data_tab <- cbind.data.frame(df4[,c('infect.at.timet','ID','HH_CERTAIN','t.index','AGE','time_risk')])
  data_t = data.table(data_tab)
  Y.df = data_t[,list(A = mean(infect.at.timet)), by = 'ID,HH_CERTAIN,t.index']
  Y.df <- setorder(Y.df,HH_CERTAIN,ID, t.index) ### To order based on t.index and not Person_ID
  
  
  #Design matrix 
  X <- df4[c('alpha0','delta0','vax2a','vax2b','time_risk','AGE1','AGE2','vax1a','vax1b','ID','HH_CERTAIN','t.index')]  
  X <- data.table(X)
  X <- setorder(X, HH_CERTAIN, ID, t.index)
  
  Y <-  Y.df$A
  
  out.list=list('Y'=Y, 'X'=X)  
  return(out.list)
}

