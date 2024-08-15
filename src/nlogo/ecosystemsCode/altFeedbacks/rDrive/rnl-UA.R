require(snowfall)
require(RNetLogo)

setwd("J:/private/SourceCode/NetLogo/fireModels/altFeedbacks/rDrive")

processors <- 7 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "logUA.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.1.0"
model.path <- "J:/private/SourceCode/NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         
n.runs <- 2000
n.params.eval <- 16
n.sv.returned <- 11

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)
params[,2] <- rep(250, n.runs)
params[,3] <- rep(10, n.runs)
params[,4] <- runif(n.runs, 0.00, 0.25)     # fire-frequency
params[,5] <- runif(n.runs, 0.1,0.6)        # frac-consumed
params[,6] <- runif(n.runs, 0, 0.95)        # seed-pred
params[,7] <- runif(n.runs, 1, 8)           # mean-ldd (patch units)        
params[,8] <- runif(n.runs, 0, 1.0)         # flamm-start
params[,9] <- runif(n.runs, 0.005, 0.15)    # base-invasion
params[,10] <- round(runif(n.runs, 1, 5))   # fire-slow
params[,11] <- runif(n.runs, 0.05, 0.20)    # fire-invasion
params[,12] <- round(runif(n.runs, 2, 8))   # base-seed-prod-3
params[,13] <- round(runif(n.runs, 2, 8))   # base-seed-prod-4
params[,14] <- round(runif(n.runs, 5, 20))   # crit-density-yng
params[,15] <- round(runif(n.runs, 5, 20))   # crit-density-old


params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- round(runif(n.runs, 0, 1))   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 1  			   # record-tag (identifier for fire record files)	
params[,18] <- rep(1.01, n.runs)       # max forest before stop
params[,19] <- as.factor(1 == 1)  # no need to write the fire record 


names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")


#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
# sfExport("ret")

r <- 1:n.runs
sfClusterEval(ls())

#print(date())
result.par <- sfClusterApplyLB(r, simfun.ua) 
#print(date())

#print(date())
#result.par <- sfClusterApplyLB(r, simfun)  # load balance?
#print(date())


# print(data.frame(t(result.par)))
op <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
op <- cbind(params, op)
names(op) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
               "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
               "crit.density.yng", "crit.density.old", "invasion", "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", 
               "n.fire", "max.fire", "sum.fire", "beta.fire")

save.image(file = "UA-v2.RData")
# system('mergeFireRecords.py')    # this binds all the firerecords together.
invisible(sfLapply(1:processors, postpro))
sfStop()
