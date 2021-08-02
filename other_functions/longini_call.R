call.longini <- function(ds.call){
##############################################################
#Model Fitting
##############################################################
inits1=list(".RNG.seed"=c(1234), ".RNG.name"='base::Wichmann-Hill')
inits2=list(".RNG.seed"=c(4567), ".RNG.name"='base::Wichmann-Hill')
inits3=list(".RNG.seed"=c(6789), ".RNG.name"='base::Wichmann-Hill')
##############################################
#Model Organization
##############################################
model_spec<-textConnection(longini_jags)
model_jags_covar<-jags.model(model_spec, 
                             inits=list(inits1,inits2, inits3),
                             data=list( 
                               'ds'=ds.call,
                               'max.hh.size'=dim(ds.call)[2],
                               'n.hh'= apply(ds.call,c(2),sum),
                               'choose_kj_mat'=choose_kj_mat
                             ),
                             n.adapt=10000, 
                             n.chains=3)

params<-c('B','Q')
##############################################
#Posterior Sampling
##############################################
posterior_samples_covar<-coda.samples(model_jags_covar, 
                                      params, 
                                      n.iter=10000)

posterior_samples.covar.all<-do.call(rbind,posterior_samples_covar)
post_median<-apply(posterior_samples.covar.all, 2, median)
sample.labs<-names(post_median)
ci<-t(hdi(posterior_samples.covar.all, credMass = 0.95))
ci <- cbind(post_median,ci)
row.names(ci)<-sample.labs
names(post_median)<-sample.labs
B.index <- grep('B',sample.labs, fixed=T)
Q.index <- grep('Q',sample.labs, fixed=T)

CPI1.mcmc <- 1- ci[B.index,]
SAR1.mcmc <- 1-ci[Q.index,]

out.list <- list('CPI1.mcmc'=CPI1.mcmc,'SAR1.mcmc')

}