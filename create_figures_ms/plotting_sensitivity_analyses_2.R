plot_VE <- function(v1,v2,v3,v1_post,v2_post,v3_post,x1,x2){
  p <- list()
  
  p[[1]] <- ggplot(v1, aes(x=factor(index), y=estimate)) +
    geom_point(color=c("black","blue","black","black"))+
    labs(tag = "A") +
    theme(plot.tag.position = c(5, 0.01))+
    # geom_text(x=0.05, y=100, label="A") +
    geom_errorbar(aes(ymin=low, ymax=up), width=0.0,color=c("black","blue","black","black"))+
    coord_cartesian(ylim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +scale_x_discrete(labels=c("1" = "0d", "2" = "1.5d","3"="3d","4"="4.5d"))+
    theme(axis.line.x=element_line(color="black"))+
    ggtitle("Partially vaccinated \n pre Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[2]] <- ggplot(v2, aes(x=factor(index), y=estimate)) +
    geom_point(color=c("black","blue","black","black"))+
    labs(tag = "B") +
    theme(plot.tag.position = c(1, 0.02))+
    geom_errorbar(aes(ymin=low, ymax=up), width=0.0,color=c("black","blue","black","black"))+
    coord_cartesian(ylim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +scale_x_discrete(labels=c("1" = "0d", "2" = "1.5d","3"="3d","4"="4.5d"))+
    theme(axis.line.x=element_line(color="black"))+
    ggtitle("Fully vaccinated \n (10-<90d after dose2) pre Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[3]] <- ggplot(v3, aes(x=factor(index), y=estimate)) +
    geom_point(color=c("black","blue","black","black"))+
    labs(tag = "C") +
    theme(plot.tag.position = c(1, 0.02))+
    geom_errorbar(aes(ymin=low, ymax=up), width=0.0,color=c("black","blue","black","black"))+
    coord_cartesian(ylim=c(-350,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +scale_x_discrete(labels=c("1" = "0d", "2" = "1.5d","3"="3d","4"="4.5d"))+
    theme(axis.line.x=element_line(color="black"))+
    ggtitle("Fully vaccinated \n (>=90d after dose2)  pre Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[4]] <- ggplot(v1_post, aes(x=factor(index), y=estimate)) +
    geom_point(color=c("black","blue","black","black"))+
    labs(tag = "D") +
    theme(plot.tag.position = c(1, 0.02))+
    geom_errorbar(aes(ymin=low, ymax=up), width=0.0,color=c("black","blue","black","black"))+
    coord_cartesian(ylim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +scale_x_discrete(labels=c("1" = "0d", "2" = "1.5d","3"="3d","4"="4.5d"))+
    theme(axis.line.x=element_line(color="black"))+
    ggtitle("Partially vaccinated \n post Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[5]] <- ggplot(v2_post, aes(x=factor(index), y=estimate)) +
    geom_point(color=c("black","blue","black","black"))+
    labs(tag = "E") +
    theme(plot.tag.position = c(1, 0.02))+
    geom_errorbar(aes(ymin=low, ymax=up), width=0.0,color=c("black","blue","black","black"))+
    coord_cartesian(ylim=c(-350,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +scale_x_discrete(labels=c("1" = "0d", "2" = "1.5d","3"="3d","4"="4.5d"))+
    theme(axis.line.x=element_line(color="black"))+
    ggtitle("Fully vaccinated \n (10-<90d after dose2) post Delta")+
    theme(plot.title = element_text(size = 12))
  
  p[[6]] <- ggplot(v3_post, aes(x=factor(index), y=estimate)) +
    geom_point(color=c("black","blue","black","black"))+
    labs(tag = "F") +
    theme(plot.tag.position = c(1, 0.02))+
    geom_errorbar(aes(ymin=low, ymax=up), width=0.0,color=c("black","blue","black","black"))+
    coord_cartesian(ylim=c(x1,x2))+
    theme_classic() +
    xlab('') +
    ylab('') +scale_x_discrete(labels=c("1" = "0d", "2" = "1.5d","3"="3d","4"="4.5d"))+
    theme(axis.line.x=element_line(color="black"))+
    ggtitle("Fully vaccinated \n (>=90d after dose2) post Delta")+
    theme(plot.title = element_text(size = 12))
  
  
  return(p)
}