# the initialization function
prepro <- function(dummy, gui, nl.path, model.path, m) {
  library(RNetLogo)
  NLStart(nl.path, nl.version=5, gui=gui)
  NLLoadModel(model.path)
}


# the simulation function
simfun.ua <- function(x) {
  
  # flush.console()
  if (x == 1) {cat(paste("Started rep 1 at: ", Sys.time(), '\n'))}
  cat(paste("Run: ", x, '\n'))
  
  NLCommand("set from-r? ", params[x, "from-r"])
  NLCommand("set run-id ", x)
  NLCommand("set burn-in-regen ", params[x, "burn-in-regen"])
  NLCommand("set max-ticks", params[x, "max-ticks"])

  NLCommand("set fire-frequency", params[x, "fire-frequency"])
  NLCommand("set fraction-consumed", params[x, "frac-consumed"])
  NLCommand("set seed-pred", params[x, "seed-pred"])
  NLCommand("set mean-ldd", params[x, "mean-ldd"])
  NLCommand("set flamm-start", params[x, "flamm-start"])
  NLCommand("set base-invasion", params[x, "base-invasion"])
  NLCommand("set fire-slow", params[x, "fire-slow"])
  NLCommand("set fire-invasion", params[x, "fire-invasion"])
  NLCommand("set base-seed-prod-3", params[x, "base-seed-prod-3"])
  NLCommand("set base-seed-prod-4", params[x, "base-seed-prod-4"])
  NLCommand("set crit-density-yng", params[x, "crit-density-yng"])
  NLCommand("set crit-density-old", params[x, "crit-density-old"])
  NLCommand("set invasion?", params[x, "invasion"])
  NLCommand("set record-tag", params[x, "record-tag"])
  NLCommand("set max-forest", params[x, "max-forest"])
  NLCommand("set write-record?", params[x, "write-record"])
  
  
  NLCommand("setup")
  
  NLCommand("go")  
  cat(paste("Year: ", NLReport("ticks"), '\n'))
  NLCommand("write-record")
  
  ret <- data.frame(c(x, NLReport("report-evaluation-variables")))
  #ret <- data.frame(x, NLReport("max-ticks"), NLReport("burn-in-regen"), NLReport("fire-frequency"), NLReport("item 3 abundances")) 
}

postpro <- function(x) {
 NLQuit()
}

