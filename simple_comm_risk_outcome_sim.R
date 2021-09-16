simple.comm.risk.outcome <- function(input_df,i,cl1){
  max.time.hh <- as.numeric(round(difftime(enddate,startdate,units="days")))+1
  input_df$max.time.hh <- max.time.hh
  df2 <-input_df
  ### Rename all of the variables on df2 as HH_CERTAIN_b, Person_ID_b, Vax_dose1b
  colnames(df2) <- paste(colnames(df2), 'b', sep='_')
  
  ### Merge df1 and df2 by HH_CERTAIN; name this df3 
  df3 <- merge(input_df,df2,by.y='HH_CERTAIN_b' , by.x='HH_CERTAIN')
  df3.spl <- split(df3, paste0(df3$ID, df3$ID_b))
  df4.spl <- lapply(df3.spl, function(x){
    x$expand.t <- NA
    x$expand.t <- x$max.time.hh#censor uninfected person at max.time.hh
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
    df.t$keep <-  T*ifelse(df.t$ID == df.t$ID_b,1,0)  
    dft.k <- df.t[df.t$keep==1,]
    dft.k <- dft.k[,c('ID','ID_b','HH_CERTAIN','vax1dose_date','vax2dose_date','vax1dose_date_b','vax2dose_date_b',
                      'infected','AGE','infected_b','max.time.hh','t.index')]
    return(dft.k)
  })
  
  df4 <-bind_rows(df4.spl)
  df4$date.current <- start+df4$t.index-1
  
  #Create elements for design matrix
  df4$delta0 <- 0
  df4$delta0[(df4$ID == df4$ID_b)] <- 1 #exogenous/community risk 
  df4$alpha0 <- 0
  df4$alpha0[df4$ID != df4$ID_b ] <- 1
  df4$vax2a <- (df4$date.current > (df4$vax2dose_date  +10))*1
  df4$vax2a[is.na(df4$vax2dose_date)] <- 0
  df4$vax2b <- (df4$date.current > (df4$vax2dose_date_b  +10))*1*df4$alpha0
  df4$vax2b[is.na(df4$vax2dose_date_b)] <- 0  
  df4$AGE1 <- 0
  df4$AGE1[(df4$AGE<60)&(df4$AGE>10)]<- 1
  df4$AGE2 <- 0
  df4$AGE2[(df4$AGE>=60)]<- 1
  df4 <- df4[!is.na(df4$HH_CERTAIN),]
  df4$time_risk <- as.numeric(pos_test_complete[df4$t.index])*df4$delta0
  #Design matrix 
  X <- df4[c('alpha0','delta0','vax2a','vax2b','time_risk','AGE1','AGE2','ID','HH_CERTAIN','t.index','date.current')]  
  X <- data.table(X)
  X <- setorder(X, HH_CERTAIN, ID,t.index)
  mat.inf <- X
  mat.inf$infected <- 0
  for (t in 1:max.time.hh){
    infect.status <-0
    log_p <- as.vector(as.matrix(X[(X$t.index==t),c('alpha0','delta0','vax2a','vax2b')]) %*% params_true) 
    ### Go back to p (probability  of transmission) with exp: 
    q <- 1 - exp(log_p)
    ##Pi needs to be a single value by ID/hhID/time point; Y should be same length
    data_tab <- cbind.data.frame(log.q=log(q),X[(X$t.index==t),])
    data_t = data.table(data_tab)
    ans = data_t[,list(A = sum(log.q)), by = 'ID,t.index']
    pi= 1- exp(ans$A)
    infect.status <- rbinom(n=length(pi),size=1,prob = pi)
    mat.inf$infected[(mat.inf$t.index==t)&(mat.inf$infected==0)]<-infect.status
  }
  df.spl <- split(mat.inf, mat.inf$ID)
  df1.spl <- lapply(df.spl, function(x){
    t_star <-ifelse(max(x$t.index[x$infected==0])==max.time.hh,1,min(x$t.index[x$infected==1])) 
    x <- x[(x$t.index ==t_star)]
    #t_noninfect <- min(x$t.index[(x$infected==0)])
    x$PCR_DATE[x$infected==0]  <- as.character(x$date.current)
    x$PCR_DATE[x$infected==1] <- as.character(x$date.current+rgamma(1,shape=2.34,scale=(1/(2.59))))#
    return(x)
  })
  df1 <-bind_rows(df1.spl)
  return(df1)
}
