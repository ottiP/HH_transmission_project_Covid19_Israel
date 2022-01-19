######### Example plotting for our sim results
plot_VE <- function(v1,v2,v3,v1_post,v2_post,v3_post,x1,x2){
  p <- list()
  
  p[[1]] <- ggplot(v1$Summary.Estimates, aes(x=ve.mean, y=index)) +
    geom_point() +
    labs(tag = "A") +
    theme(plot.tag.position = c(5, 0.01))+
    geom_errorbar(aes(xmin=ve.lcl, xmax=ve.ucl), width=0.0)+
    coord_cartesian(xlim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +
    geom_vline(xintercept=0, lty=2, col='gray') +
    geom_point(aes(x=v1$Individual.VE.Estimates[1], y=0), col='red', cex=3 )+
    geom_errorbar(aes(y=0, xmin=v1$Individual.VE.Estimates[2], xmax=v1$Individual.VE.Estimates[3]), width=0.0, col='red', lwd=2) + 
    ggtitle("Partially vaccinated \n pre Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[2]] <- ggplot(v2$Summary.Estimates, aes(x=ve.mean, y=index)) +
    geom_point() +
    geom_errorbar(aes(xmin=ve.lcl, xmax=ve.ucl), width=0.0)+
    labs(tag = "B") +
    theme(plot.tag.position = c(5, 0.01))+
    coord_cartesian(xlim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +
    geom_vline(xintercept=0, lty=2, col='gray') +
    geom_point(aes(x=v2$Individual.VE.Estimates[1], y=0), col='red', cex=3 )+
    geom_errorbar(aes(y=0, xmin=v2$Individual.VE.Estimates[2], xmax=v2$Individual.VE.Estimates[3]), width=0.0, col='red', lwd=2)+
    ggtitle("Fully vaccinated \n (10-<90d after dose2) pre Delta")+
    theme(plot.title = element_text(size = 12))
  
  
  p[[3]] <- ggplot(v3$Summary.Estimates, aes(x=ve.mean, y=index)) +
    geom_point() +
    labs(tag = "C") +
    theme(plot.tag.position = c(5, 0.01))+
    geom_errorbar(aes(xmin=ve.lcl, xmax=ve.ucl), width=0.0)+
    coord_cartesian(xlim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +
    geom_vline(xintercept=0, lty=2, col='gray') +
    geom_point(aes(x=v3$Individual.VE.Estimates[1], y=0), col='red', cex=3 )+
    geom_errorbar(aes(y=0, xmin=v3$Individual.VE.Estimates[2], xmax=v3$Individual.VE.Estimates[3]), width=0.0, col='red', lwd=2)+
    ggtitle("Fully vaccinated \n (>=90d after dose2)  pre Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[4]] <- ggplot(v1_post$Summary.Estimates, aes(x=ve.mean, y=index)) +
    geom_point() +
    geom_errorbar(aes(xmin=ve.lcl, xmax=ve.ucl), width=0.0)+
    coord_cartesian(xlim=c(x1,x2))+
    labs(tag = "D") +
    theme(plot.tag.position = c(5, 0.01))+
    theme_classic() +
    xlab('') +
    ylab('') +
    geom_vline(xintercept=0, lty=2, col='gray') +
    geom_point(aes(x=v1_post$Individual.VE.Estimates[1], y=0), col='red', cex=3 )+
    geom_errorbar(aes(y=0, xmin=v1_post$Individual.VE.Estimates[2], xmax=v1_post$Individual.VE.Estimates[3]), width=0.0, col='red', lwd=2)+
    ggtitle("Partially vaccinated \n post Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[5]] <- ggplot(v2_post$Summary.Estimates, aes(x=ve.mean, y=index)) +
    geom_point() +
    geom_errorbar(aes(xmin=ve.lcl, xmax=ve.ucl), width=0.0)+
    labs(tag = "E") +
    theme(plot.tag.position = c(5, 0.01))+
    coord_cartesian(xlim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +
    geom_vline(xintercept=0, lty=2, col='gray') +
    geom_point(aes(x=v2_post$Individual.VE.Estimates[1], y=0), col='red', cex=3 )+
    geom_errorbar(aes(y=0, xmin=v2_post$Individual.VE.Estimates[2], xmax=v2_post$Individual.VE.Estimates[3]), width=0.0, col='red', lwd=2)+
    ggtitle("Fully vaccinated \n (10-<90d after dose2) post Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[6]] <- ggplot(v3_post$Summary.Estimates, aes(x=ve.mean, y=index)) +
    geom_point() +
    geom_errorbar(aes(xmin=ve.lcl, xmax=ve.ucl), width=0.0)+
    labs(tag = "F") +
    theme(plot.tag.position = c(5, 0.01))+
    coord_cartesian(xlim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +
    geom_vline(xintercept=0, lty=2, col='gray') +
    geom_point(aes(x=v3_post$Individual.VE.Estimates[1], y=0), col='red', cex=3 )+
    geom_errorbar(aes(y=0, xmin=v3_post$Individual.VE.Estimates[2], xmax=v3_post$Individual.VE.Estimates[3]), width=0.0, col='red', lwd=2)+
    ggtitle("Fully vaccinated \n (>=90d after dose2) post Delta")+
    theme(plot.title = element_text(size = 12))
  
  return(p)
}
