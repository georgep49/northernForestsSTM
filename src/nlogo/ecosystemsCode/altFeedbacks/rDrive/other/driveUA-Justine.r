library(RNetlogo)
setwd("C:/Users/gper020/Desktop/Trade-off analysis")
m.path <- "C:/Users/gper020/Desktop/Trade-off analysis/Final model version 1.3.nlogo"
nl.path <- "C:/Program Files (x86)/NetLogo 5.0.4"

source(file = "driversRNL-Stick.r")
source (file = "C:/Users/gper020/Desktop/Trade-off analysis/brokenStickFunction.r")

# Load the model (if you don't want to see it, set show = FALSE)
init.rnl(nl.path, m.path, show = TRUE)

#remove hypercube analysis
# library(lhs)
# n <- 10000      # number of simulations
# k <- 4         # number of parameters being assessed

# # Build the hypercube here
# L <- randomLHS(n = n, k = k)
# 
# # This rescales the parameters of interest
# speed.p <- as.integer(qunif(L[,1], 1, 6)) # 1
# view.p <- as.integer(qunif(L[,2], 1, 6)) # 5
# mortality.p <- qunif(L[,3], 0.0001, 0.01)#many
# foraging.p <-qunif(L[,4],0, 0.5)
# 
# # This is a matrix of the parameter values to use
# p.mtx <- matrix(c(speed.p, view.p, mortality.p, foraging.p), ncol = 4)
# 
# # Names in the matrix MUST match variable names in NetLogo
# colnames(p.mtx) <- c("speed", "view-distance", "per-step-mortality","p-start")

#source brokenStick matrix

p.mtx <- allocate.points(reps = 10, n.traits = 5, total.points = 100 )$stick.mtx

# Do the analysis
model.ua3 <- sim.ua.stick(p.mtx, save.to.file = TRUE, write.file = "LHS output3.csv")

# If you wanted to split
#model.ua1 <- sim.ua(p.mtx[1:100,])
#model.ua2 <- sim.ua(p.mtx[101:200,])


                     
                     
                     
                     
                     