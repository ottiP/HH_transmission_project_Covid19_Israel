compute_CI_comb <- function(Coefs, Vcov, parmnames){
  colnames(Vcov) <- names(Coefs)
  rownames(Vcov) <- names(Coefs)
  comb.effect <- sum(Coefs[parmnames])
  Var <- diag(Vcov)
  names(Var) <- colnames(Vcov)
  if(length(parmnames)>1){
    parm.combos <- combn(parmnames,2,simplify=F)
    cov.piece <- sapply(parm.combos, function(zz){
      Cov <- 2 * Vcov[zz[1], zz[2] ]
    })
    var.comb.effect <- sum(Var[parmnames]) + sum(cov.piece)
  }else{
    var.comb.effect <- Var[parmnames]
  }
  lcl <- comb.effect - 1.96*sqrt(var.comb.effect)
  ucl <- comb.effect + 1.96*sqrt(var.comb.effect)
  parm.combos.ci <- c('mean'=comb.effect, 'lcl'=lcl,'ucl'=ucl)
  out.list <- list('parm.combos.ci'=parm.combos.ci, 'var.comb.effect'=var.comb.effect)
  return(out.list)
}


compute_law_tot_var <- function(RR,SE){
  SE_var <- sqrt(var(unlist(RR))+mean(SE*SE))
  return(SE_var)
}


combined.estimates <- function(Set.Coefs=parameters,Set.Vcov=OI,Set.parmnames=c('vax1_contact2')){
  comb1 <- mapply(FUN=compute_CI_comb, Coefs=Set.Coefs, Vcov=Set.Vcov, MoreArgs=list(parmnames=Set.parmnames), SIMPLIFY=F )
  var.comb1 <- sapply(comb1, '[[', 'var.comb.effect')
  est.comb1 <- as.data.frame(t(sapply(comb1, '[[', 'parm.combos.ci')))
  est.comb1$index <- 1:nrow(est.comb1)
  tot.se <-  sqrt(var(est.comb1[,'mean']) + mean(var.comb1))
  
  combine.mean <- mean(est.comb1[,'mean'])
  combine.lci <- combine.mean - 1.96*tot.se
  combine.uci <- combine.mean + 1.96*tot.se
  
  combine.ci <- c(combine.mean,combine.lci, combine.uci)
  combine.ve.ci <- 100*(1-exp(combine.ci))
  combine.rr.ci <- exp(combine.ci)
  
  est.comb1$ve.mean <- 100*(1-exp(est.comb1[,1]))
  est.comb1$ve.lcl <- 100*(1-exp(est.comb1[,2]))
  est.comb1$ve.ucl <- 100*(1-exp(est.comb1[,3]))
  
  p1 <- ggplot(est.comb1, aes(x=ve.mean, y=index)) +
    geom_point() +
    geom_errorbar(aes(xmin=ve.lcl, xmax=ve.ucl), width=0.0)+
    coord_cartesian(xlim=c(-100,100))+
    theme_classic() +
    xlab('Vaccine Effectiveness (%)') +
    ylab('Draw from delay distribution') +
    geom_vline(xintercept=0, lty=2, col='gray') +
    geom_point(aes(x=combine.ve.ci[1], y=0), col='red', cex=3 )+
    geom_errorbar(aes(y=0, xmin=combine.ve.ci[2], xmax=combine.ve.ci[3]), width=0.0, col='red', lwd=2)
  
  out.list <- list('catterplot'=p1, 'Summary.Estimates'=est.comb1, 'Individual.Estimates'=combine.ci,'Individual.VE.Estimates'=combine.ve.ci )
}
