library(nlrx)
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.4.0")
modelpath <- file.path("src/nlogo/model.nlogo")
outpath <- file.path("E:/dump")

# Setup nl object
nl <- nl(nlversion = "6.4.0",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

# Attach experiment
nl@experiment <- experiment(expname="wolf-sheep",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            runtime=1,
                            evalticks=1,
                            metrics=c("count turtles"),
                            variables = list('x' = list(min=50, max=150, qfun="qunif")))

# Attach simdesign
nl@simdesign <- simdesign_lhs(nl=nl,
                               samples=2,
                               nseeds=3,
                               precision=3)

# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

# Run all simulations (loop over all siminputrows and simseeds)
results <- run_nl_all(nl)

# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Write output to outpath of experiment within nl
write_simoutput(nl)

# Do further analysis:
analyze_nl(nl)