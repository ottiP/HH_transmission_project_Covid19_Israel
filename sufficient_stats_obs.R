

a1 <- readRDS(sim.data.df,file=paste("./Data/simulated_data.rds"))

a1$agec <- NA
a1$agec[a1$AGE>=0 & a1$AGE<40] <- 1
a1$agec[a1$AGE>=40 & a1$AGE<65] <- 2
a1$agec[a1$AGE>=65 & a1$AGE<110] <- 3

