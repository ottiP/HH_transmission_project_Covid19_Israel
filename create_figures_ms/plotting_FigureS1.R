d2 <- readRDS('A:\\intermediate_data\\clean.hh.data.rds')
d2<-d2[!(d2$HH_CERTAIN=='NULL'),]
pos_date <- d2 %>%
  group_by(PCR_DATE) %>%
  summarise(pos_test = sum(IS_POSITIVE_CD==2))
pos_date<-pos_date[109:531,]
cases<- as.data.frame(pos_date)
cases$dates <-as.Date(cases$PCR_DATE)
cases$cases_av <-mov_av_data$pos_test
cases$cases_av_st<-scale(log(mov_av_data$pos_test))
g <- list()
g[[1]]<-ggplot(cases,aes(x=dates,y=cases$pos_test),group=1)+geom_line()+ylab("Positive PCR test")+xlab("")+geom_vline(xintercept=as.numeric(as.Date("2020-12-20")), linetype=3,color='red')+geom_vline(xintercept=as.numeric(as.Date("2021-06-01")), linetype=3)+geom_line(aes(x=cases$dates,y=cases$cases_av),color='blue')+theme_bw()
# +geom_line(aes(x=cases$dates,y=cases$cases_av_st),color='orange')
g[[2]]<-ggplot(cases,aes(x=dates,y=cases$cases_av_st),group=1)+geom_line()+ylab("Standardized Positive PCR test")+xlab("")+geom_vline(xintercept=as.numeric(as.Date("2020-12-20")), linetype=3,color='red')+geom_vline(xintercept=as.numeric(as.Date("2021-06-01")), linetype=3)+theme_bw()
grid.arrange(grobs=g,ncol=2,nrow=1,fontsize=18,heights=c(4), widths=c(1,1))
ggsave("Figure_S1.pdf",dpi=300,arrangeGrob(grobs=g,ncol=2,nrow=1,width=0.5, height=10))
