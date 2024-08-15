#####################################
# Scenario 1:: Seed predation rates
#####################################

setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
setwd('/home/gper020/Scifac/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 400
n.reps <- 20

n.params.eval <- 19
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)				# called from R?
params[,2] <- rep(5000, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length
params[,4] <- rep(0.00, n.runs)     		# fire-frequency
params[,5] <- rep(0.4, n.runs)        		# frac-consumed
params[,6] <- rep(seq(0, 0.95, 0.05), times = n.reps)        # seed-pred
params[,7] <- rep(4, n.runs)            # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)          # flamm-start
params[,9] <- rep(0.05, n.runs)          # base-invasion
params[,10] <- round(rep(2, n.runs))    # fire-slow
params[,11] <- rep(0.10, n.runs)        # fire-invasion
params[,12] <- round(rep(4, n.runs))    # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))    # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old
		

params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- round(0)   # invasion?
params[,16] <- as.factor(i.test == 1)  # invasion?
params[,17] <- 1				   # record-tag (identifier for fire record files)	
params[,18] <- rep(0.95, n.runs)       # max forest before stop
params[,19] <- as.factor(i.test == 1)  # no need to write the fire record 


names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves

cat("Into scenario 1...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
pred.scenario <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
pred.scenario <- cbind(params, pred.scenario)
names(pred.scenario) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
               "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
               "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
               "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
               "sum.fire", "ave.fire")

save.image(file = "predationScenario.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()



#####################################
# Scenario 2:: Consumption rates
#####################################
rm(list = ls())

setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log2.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 420
n.reps <- 20

n.params.eval <- 19
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)				# called from R?
params[,2] <- rep(5000, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length
params[,4] <- rep(0.00, n.runs)     		# fire-frequency
params[,5] <- rep(seq(0, 1.0, 0.05), times = n.reps)        		# frac-consumed
params[,6] <- rep(0.00, n.runs)        # seed-pred
params[,7] <- rep(4, n.runs)           # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)         # flamm-start
params[,9] <- rep(0.05, n.runs)    # base-invasion
params[,10] <- round(rep(2, n.runs))   # fire-slow
params[,11] <- rep(0.10, n.runs)    # fire-invasion
params[,12] <- round(rep(4, n.runs))   # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))   # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old
		

params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- round(0)   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 2				# record-tag (identifier for fire record files)	
params[,18] <- rep(0.95, n.runs)    # max forest before stop
params[,19] <- as.factor(i.test == 1)  # no need to write the fire record 

n.params.eval <- ncol(params)

names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
cat("Into scenario 2...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
cons.scenario <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
cons.scenario <- cbind(params, cons.scenario)
names(cons.scenario) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                          "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                          "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                          "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                          "sum.fire", "ave.fire")

save.image(file = "consumptionScenario.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()



#####################################
# Scenario 3:: consumption x LDD rates
#####################################
rm(list = ls())

setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log3.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 2200
n.reps <- 20

n.params.eval <- 19
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)				# called from R?
params[,2] <- rep(5000, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length
params[,4] <- rep(0.0, n.runs)    		# fire-frequency

cmbn <- expand.grid(seq(0, 1.0, 0.1), seq(0, 0.9, 0.1)) 
params[,5] <- rep(as.matrix(cmbn[1]), times = n.reps)        	# frac-consumed
params[,6] <- rep(as.matrix(cmbn[2]), times = n.reps)            # seed-pred

params[,7] <- rep(4, n.runs)           # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)         # flamm-start
params[,9] <- rep(0.05, n.runs)    # base-invasion
params[,10] <- round(rep(2, n.runs))   # fire-slow
params[,11] <- rep(0.10, n.runs)    # fire-invasion
params[,12] <- round(rep(4, n.runs))   # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))   # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old
		

params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- round(0)   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 3				# record-tag (identifier for fire record files)	
params[,18] <- rep(0.95, n.runs)    # max forest before stop
params[,19] <- as.factor(i.test == 1)  # no need to write the fire record 

n.params.eval <- ncol(params)

names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
cat("Into scenario 3...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
conspred.scenario <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
conspred.scenario  <- cbind(params, conspred.scenario )
names(conspred.scenario ) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                               "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                               "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                               "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                               "sum.fire", "ave.fire")

save.image(file = "consumptionByPredScenario.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()






#####################################
# Scenario 4:: fire frequency
#####################################
rm(list = ls())

setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log4.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 1040
n.reps <- 20

n.params.eval <- 19
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)				# called from R?
params[,2] <- rep(500, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length
params[,4] <- rep(seq(0, 0.25, 0.01), times = n.reps)    		# fire-frequency
params[,5] <- rep(0.4, n.runs)        	# frac-consumed
params[,6] <- rep(0.0, n.runs)        # seed-pred
params[,7] <- rep(4, n.runs)           # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)         # flamm-start
params[,9] <- rep(0.05, n.runs)    # base-invasion
params[,10] <- round(rep(2, n.runs))   # fire-slow
params[,11] <- rep(0.10, n.runs)    # fire-invasion
params[,12] <- round(rep(4, n.runs))   # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))   # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old
		

params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- rep(each = 520, c(0,1))   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 4				# record-tag (identifier for fire record files)	
params[,18] <- rep(1.01, n.runs)    # max forest before stop
params[,19] <- as.factor(1 == 1)  #  need to write the fire record 

n.params.eval <- ncol(params)

names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
cat("Into scenario 4...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
ffreq.scenario  <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
ffreq.scenario <- cbind(params, ffreq.scenario)
names(ffreq.scenario) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                           "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                           "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                           "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                           "sum.fire", "ave.fire")

save.image(file = "fireFreq.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()


#####################################
# Scenario 5:: fire - invasion - rats
#####################################
rm(list = ls())

setwd('~/Scifac/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log5e.txt")
sfClusterSetupRNG()

gui <- FALSE

#nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
#model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"

nl.path <- "/home/gper020/bin/netlogo-5.1.0"
model.path <- "/home/gper020/Scifac/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"


#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 5200
n.reps <- 20

n.params.eval <- 18
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)				# called from R?
params[,2] <- rep(500, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length

cmbn <- expand.grid(seq(0.0, 0.25, 0.01), 0.8, c(0.0, 0.2, 0.4, 0.6, 0.8)) 
params[,4] <- rep(as.matrix(cmbn[1]), times = n.reps * 2)    		# fire-frequency [ * 2 for invasion]
params[,5] <- rep(as.matrix(cmbn[2]), times = n.reps * 2)       # frac-consumed [ * 2 for invasion]
params[,6] <- rep(as.matrix(cmbn[3]), times = n.reps * 2)       # seed-pred [ * 2 for invasion]

params[,7] <- rep(4, n.runs)           # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)         # flamm-start
params[,9] <- rep(0.05, n.runs)    # base-invasion
params[,10] <- round(rep(2, n.runs))   # fire-slow
params[,11] <- rep(0.10, n.runs)    # fire-invasion
params[,12] <- round(rep(4, n.runs))   # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))   # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old
		

params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- rep(each = 2600, c(0,1))   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 5				# record-tag (identifier for fire record files)	
params[,18] <- rep(1.01, n.runs)    # max forest before stop
params[,19] <- as.factor(1 == 1)  #  need to write the fire record 

n.params.eval <- ncol(params)

names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
cat("Into scenario 5...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
fir5e.scenario <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
fir5e.scenario <- cbind(params, fir5e.scenario)
names(fir5e.scenario) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                         "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                         "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                         "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                         "sum.fire", "ave.fire")

save.image(file = "FIR5e.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()

#####################################
# Scenario 5b:: fire - invasion - rats (bolster 5 with seed pred = 0.2, 0.8)
#####################################
rm(list = ls())

setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log5b.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 4160
n.reps <- 20

n.params.eval <- 18
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)  			# called from R?
params[,2] <- rep(500, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length

cmbn <- expand.grid(seq(0, 0.25, 0.01), c(0.2, 0.4), c(0.2, 0.8)) 
params[,4] <- rep(as.matrix(cmbn[1]), times = n.reps * 2)    		# fire-frequency [ * 2 for invasion]
params[,5] <- rep(as.matrix(cmbn[2]), times = n.reps * 2)        	# frac-consumed [ * 2 for invasion]
params[,6] <- rep(as.matrix(cmbn[3]), times = n.reps * 2)            # seed-pred [ * 2 for invasion]

params[,7] <- rep(4, n.runs)           # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)         # flamm-start
params[,9] <- rep(0.05, n.runs)    # base-invasion
params[,10] <- round(rep(2, n.runs))   # fire-slow
params[,11] <- rep(0.10, n.runs)    # fire-invasion
params[,12] <- round(rep(4, n.runs))   # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))   # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old


params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- rep(each = 2080, c(0,1))   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 6				# record-tag (identifier for fire record files)	
params[,18] <- rep(1.01, n.runs)    # max forest before stop
params[,19] <- as.factor(1 == 1)  #  need to write the fire record 

n.params.eval <- ncol(params)

names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
cat("Into scenario 5b...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
fir5b.scenario <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
fir5b.scenario <- cbind(params, fir5b.scenario)
names(fir5b.scenario) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                         "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                         "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                         "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                         "sum.fire", "ave.fire")

save.image(file = "FIR5b.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()

#####################################
# Scenario 5c:: fire - invasion - rats (bolster 5a and b with seed pred = 0.4)
#####################################
rm(list = ls())

setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log5c.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 2080
n.reps <- 20

n.params.eval <- 18
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)    		# called from R?
params[,2] <- rep(500, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length

cmbn <- expand.grid(seq(0, 0.25, 0.01), c(0.2, 0.4), 0.4) 
params[,4] <- rep(as.matrix(cmbn[1]), times = n.reps * 2)    		# fire-frequency [ * 2 for invasion]
params[,5] <- rep(as.matrix(cmbn[2]), times = n.reps * 2)        	# frac-consumed [ * 2 for invasion]
params[,6] <- rep(as.matrix(cmbn[3]), times = n.reps * 2)            # seed-pred [ * 2 for invasion]

params[,7] <- rep(4, n.runs)           # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)         # flamm-start
params[,9] <- rep(0.05, n.runs)    # base-invasion
params[,10] <- round(rep(2, n.runs))   # fire-slow
params[,11] <- rep(0.10, n.runs)    # fire-invasion
params[,12] <- round(rep(4, n.runs))   # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))   # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old


params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- rep(each = 1040, c(0,1))   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 7				# record-tag (identifier for fire record files)	
params[,18] <- rep(1.01, n.runs)    # max forest before stop
params[,19] <- as.factor(1 == 1)  #  need to write the fire record 

n.params.eval <- ncol(params)

names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
cat("Into scenario 5c...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
fir5c.scenario <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
fir5c.scenario <- cbind(params, fir5c.scenario)
names(fir5c.scenario) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                           "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                           "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                           "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                           "sum.fire", "ave.fire")

save.image(file = "FIR5c.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()

#####################################
# Scenario 5d:: fire - invasion - rats (bolster 5a to c with seed mvt = 0.6 and 0.8:: in line with NZ data)
#####################################
rm(list = ls())

setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log5d.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 10400
n.reps <- 20

n.params.eval <- 18
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)      	# called from R?
params[,2] <- rep(500, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length

cmbn <- expand.grid(seq(0, 0.25, 0.01), c(0.6, 0.8), c(0, 0.2, 0.4, 0.6, 0.8)) 
params[,4] <- rep(as.matrix(cmbn[1]), times = n.reps * 2)    		# fire-frequency [ * 2 for invasion]
params[,5] <- rep(as.matrix(cmbn[2]), times = n.reps * 2)       # frac-consumed [ * 2 for invasion]
params[,6] <- rep(as.matrix(cmbn[3]), times = n.reps * 2)       # seed-pred [ * 2 for invasion]

params[,7] <- rep(4, n.runs)           # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)         # flamm-start
params[,9] <- rep(0.05, n.runs)    # base-invasion
params[,10] <- round(rep(2, n.runs))   # fire-slow
params[,11] <- rep(0.10, n.runs)    # fire-invasion
params[,12] <- round(rep(4, n.runs))   # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))   # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old


params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- rep(each = 2600, c(0,1))   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 7				# record-tag (identifier for fire record files)	
params[,18] <- rep(1.01, n.runs)    # max forest before stop
params[,19] <- as.factor(1 == 1)  #  need to write the fire record 

n.params.eval <- ncol(params)

names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
cat("Into scenario 5d...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
fir5d.scenario <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
fir5d.scenario <- cbind(params, fir5d.scenario)
names(fir5d.scenario) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                           "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                           "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                           "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                           "sum.fire", "ave.fire")

save.image(file = "FIR5d.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()

#####################################
# Scenario 1b:: Seed predation rates
#####################################

#setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
# setwd('~/windows/J/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
setwd("~/Scifac/private/SourceCode/NetLogo/fireModels/altFeedbacks/rDrive")

require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "logXoom2.txt")
sfClusterSetupRNG()

gui <- FALSE

#nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
nl.path <- "/home/gper020/bin/netlogo-5.1.0"

#model.path <- "~/windows/J/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
model.path <- "/home/gper020/Scifac/private/SourceCode/NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"

#model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"

#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 600
n.reps <- 20

n.params.eval <- 19
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)  			# called from R?
params[,2] <- rep(5000, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length
params[,4] <- rep(0.00, n.runs)     		# fire-frequency
params[,5] <- rep(0.4, n.runs)        		# frac-consumed
params[,6] <- rep(seq(0.6, 0.745, 0.005), times = n.reps)        # seed-pred
params[,7] <- rep(4, n.runs)            # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)          # flamm-start
params[,9] <- rep(0.05, n.runs)          # base-invasion
params[,10] <- round(rep(2, n.runs))    # fire-slow
params[,11] <- rep(0.10, n.runs)        # fire-invasion
params[,12] <- round(rep(4, n.runs))    # base-seed-prod-3
params[,13] <- round(rep(4, n.runs))    # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old


params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- round(0)   # invasion?
params[,16] <- as.factor(i.test == 1)  # invasion?
params[,17] <- 1				   # record-tag (identifier for fire record files)	
params[,18] <- rep(0.95, n.runs)       # max forest before stop
params[,19] <- as.factor(i.test == 1)  # no need to write the fire record 


names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves

cat("Into scenario 1...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
pred.zoom.scenario2 <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
pred.zoom.scenario2 <- cbind(params, pred.zoom.scenario2)
names(pred.zoom.scenario2) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                          "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                          "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                          "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                          "sum.fire", "ave.fire")

save.image(file = "predationZoomScenario2.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()

#####################################
# Scenario 2b:: Consumption rates / low avail
#####################################
rm(list = ls())

setwd('J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/rDrive/')
require(snowfall)
require(RNetLogo)

processors <- 6 # detectCores()
sfInit(parallel=TRUE, cpus=processors, type="SOCK", slaveOutfile = "log2b.txt")
sfClusterSetupRNG()

gui <- FALSE

nl.path <- "C:/Program Files (x86)/NetLogo 5.0.5"
model.path <- "J:/private/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
#model.path <- "C:/Users/gper020/Documents/JDrive/SourceCode/Model Code NetLogo/fireModels/altFeedbacks/altFeedbacks-Regen.nlogo"
source("analyses-UA.r")

invisible(sfLapply(1:processors, prepro, gui = gui, nl.path = nl.path, model.path = model.path))         

n.runs <- 420
n.reps <- 20

n.params.eval <- 19
n.sv.returned <- 12

params <- matrix(rep(0, n.runs * n.params.eval), nc = n.params.eval)
params[,1] <- rep(TRUE, n.runs)  			# called from R?
params[,2] <- rep(5000, n.runs)				# length
params[,3] <- rep(10, n.runs)				# burn-in length
params[,4] <- rep(0.00, n.runs)     		# fire-frequency
params[,5] <- rep(seq(0, 1.0, 0.05), times = n.reps)        		# frac-consumed
params[,6] <- rep(0.00, n.runs)        # seed-pred
params[,7] <- rep(4, n.runs)           # mean-ldd (patch units)        
params[,8] <- rep(0.5, n.runs)         # flamm-start
params[,9] <- rep(0.05, n.runs)    # base-invasion
params[,10] <- round(rep(2, n.runs))   # fire-slow
params[,11] <- rep(0.10, n.runs)    # fire-invasion
params[,12] <- round(rep(2, n.runs))   # base-seed-prod-3
params[,13] <- round(rep(2, n.runs))   # base-seed-prod-4
params[,14] <- round(rep(10, n.runs))   # crit-density-yng
params[,15] <- round(rep(10, n.runs))   # crit-density-old


params <- as.data.frame(params)   # do this before booleans, etc. 

i.test <- round(0)   # invasion?
params[,16] <- as.factor(i.test == 1)               # invasion?
params[,17] <- 2				# record-tag (identifier for fire record files)	
params[,18] <- rep(0.95, n.runs)    # max forest before stop
params[,19] <- as.factor(i.test == 1)  # no need to write the fire record 

n.params.eval <- ncol(params)

names(params) <- c("from-r", "max-ticks", "burn-in-regen", "fire-frequency", "frac-consumed", "seed-pred",
                   "mean-ldd", "flamm-start", "base-invasion", "fire-slow", "fire-invasion", "base-seed-prod-3", "base-seed-prod-4",
                   "crit-density-yng", "crit-density-old", "invasion", "record-tag", "max-forest", "write-record")

#ret.params <- 6   # Number of parameters to return from Netlogo (length of report-evaluation-variables list)
#ret <- matrix(rep(-1, n.runs * ret.params), nc = ret.params)

sfExport("params")   # preload the parameters file onto the slaves
cat("Into scenario 2...")

r <- 1:n.runs
sfClusterEval(ls())
result.par <- sfClusterApplyLB(r, simfun.ua) 



# print(data.frame(t(result.par)))
consLowProd.scenario <- data.frame(matrix(unlist(t(result.par)), nc = n.sv.returned, byrow = TRUE))
consLowProd.scenario <- cbind(params, consLowProd.scenario)
names(consLowProd.scenario) <- c("from.r", "max.ticks", "burn.in.regen", "fire.frequency", "frac.consumed", "seed.pred",
                          "mean.ldd", "flamm.start", "base.invasion", "fire.slow", "fire.invasion", "base.seed.prod3", "base.seed.prod.4",
                          "crit.density.yng", "crit.density.old", "invasion", "record.tag", "max.forest", "write.record",
                          "run.id", "abund1", "abund2", "abund3", "abund4", "abund5", "bft", "ticks",  "n.fire", "max.fire", 
                          "sum.fire", "ave.fire")

save.image(file = "consumptionLowProdScenario.RData")
cat(paste(" done and saved.", '\n'))

# system('mergeFireRecords.py')    # this binds all the fire-records together.
invisible(sfLapply(1:processors, postpro))
sfStop()
