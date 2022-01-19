compute_AIC <- function(ll,p){
  
  AIC <- 2*p+2*ll
  return(AIC)
  
}

result1<- readRDS("../Results/result100_Jul_nointerA.rds")
min_ll1 <- sapply(result1,'[[','minimum')
parms1 <- sapply(result1,'[[','estimate')
p1 <- nrow(parms1)
AIC1 <- compute_AIC(min_ll1,p1)

result2<- readRDS("../Results/result100_Jul_allinter_TestA.rds")
min_ll2 <- sapply(result2,'[[','minimum')
parms2 <- sapply(result2,'[[','estimate')
p2 <- nrow(parms2)
AIC2 <- compute_AIC(min_ll2,p2)




pdf("hist_12.pdf")
hist(AIC1-AIC2, col="gray", xlab="AIC difference",ylab="",main=paste("AIC difference between  model 1 and model 2"))
dev.off()

pdf("hist_13.pdf")
hist(AIC2-AIC3, col="gray", xlab="AIC difference",ylab="",main=paste("AIC difference between  model 2 and model 3"))
dev.off()

