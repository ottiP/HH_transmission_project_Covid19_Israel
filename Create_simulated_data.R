##### Create the data set

##### Generate the synthetic data with 10k HHs and store them as a dataframe: the data is already restricted in the study 
##### period: 15th June 2020 - 24 March 2021
N.HH <- 40000
sim.data.ls <- pblapply(1:N.HH, gen.hh,CPI=(1-0.9995), prob.trans.day=(1-0.968),irr.vax1=0.5,irr.vax2=1)
#The format is the same as the real data, though the values are completely fake and relevant only to play with the code
n.infect.hh <- sapply(sim.data.ls,function(x) sum(x$infected))
sim.data.pos.hh <- sim.data.ls[(n.infect.hh>0)]
sim.data.df <- bind_rows(sim.data.pos.hh)
saveRDS(sim.data.df,file=paste("../Data/simulated_data.rds"),compress=TRUE)
