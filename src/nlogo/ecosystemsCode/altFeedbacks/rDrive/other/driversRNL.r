# This loads the model
init.rnl <- function(nl.path, mod.path, show = FALSE)
{
  library(RNetLogo)
  # nl.path <- "C:/Program Files/NetLogo 5.0.4/"
  NLStart(nl.path, gui = show, nl.version = 5)
  
  # model.path <- "C:/Users/jatk024/Desktop/Final model version 1.3.nlogo"
  NLLoadModel(mod.path)
}


# This does all the uncertainty analysis work
sim.ua <- function(param.mtx, save.to.file = FALSE, write.file = "dump.csv") 
{
  require(RNetLogo)  
  
  n.params <- dim(param.mtx)[2]
  n.sims <- dim(param.mtx)[1]	
  
  # This is the matrix that stores the outputs.  Set ncol to the number of outputs.
  # This will be the length of the list returned by report-evaluation-variables
  ret <- matrix(NA, nrow = n.sims, ncol = 4)
  
  pb <- txtProgressBar(min = 0, max = n.sims, initial = 0, style = 3)
  
  n <- 1
  while (n <= n.sims)
  {
    NLCommand("reset-to-baseline")
    
    base.speed <- 1
    base.perception <- 1
    base.mortality <- 0.01
    base.foraging <- 0
    base.habitat <- 50
    
    speed.current <- base.speed + (param.mtx[n, 1] * 0.04)
    perception.current <- base.perception + (param.mtx[n, 2] * 0.04)
    mortality.current <- base.mortality - (param.mtx[n, 3] * 0.000099)
    foraging.current <- base.foraging + (param.mtx[n, 4] * 0.005)
    habitat.current <- base.habitat - (param.mtx[n, 5] * 0.49)
    
    p <- 1
    
    while (p <= n.params)
    {
      cmd <- paste("set", colnames(param.mtx)[p], param.mtx[n,p], sep = " ")
      NLCommand(cmd)
      p <- p + 1
    }
    
    setTxtProgressBar(pb, n)
    NLCommand("setup")
    NLCommand("set speed", speed.current)
    NLCommand("set view-distance", perception.current)
    NLCommand("set per-step-mortality", mortality.current)
    NLCommand("set p-start", foraging.current)
    NLCommand("set min-habitat-size", habitat.current)
    #NLDoCommandWhile("count groups > 0", "go")
    NLCommand("go")
    #NLCommand("final-tidy")
    ret[n,] <- NLReport("report-evaluation-variables")
    
    n <- n + 1
  }
  
  close(pb)
  # Change this if you change the things being sent back from NetLogo
  colnames(ret) <- c("n.patches", "success", "mean steps", "st dev steps")
  
  if (save.to.file == TRUE)
  {
      f <- paste("C:/Users/jatk024/Desktop/", write.file, sep = "") 
      write.csv(ret, file = f)
  }
  
  
  list(ua = as.data.frame(ret), inputs = as.data.frame(param.mtx))
}
