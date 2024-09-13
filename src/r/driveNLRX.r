library(nlrx)

# This needs NetLogo 6.3.0 or earlier... to do with paths to extensions
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.3.0")
modelpath <- file.path("H://files.fos/private/Research/NgaRakauTaketake/stm/src/nlogo/stmNorthernForests_63.nlogo")
outpath <- file.path("C:/Temp")


# Setup nl object
nl <- nl(nlversion = "6.3.0",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

# Attach experiment
nl@experiment <- experiment(expname="test",
                            idrunnum = "nlrx-info",
                            outpath=outpath,
                            repetition=1,
                            tickmetrics="true",
                            idsetup="setup",
                            idgo="go",
                            runtime=50,
                            evalticks=seq(40,50),
                            metrics=c("abundances"),
                            variables = list(),
                            constants = list("fire-frequency" = 0.15))

# Attach simdesign
nl@simdesign <- simdesign_simple(nl = nl, nseeds = 3)

# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

# Run all simulations (loop over all siminputrows and simseeds)
#library(future)
# plan(multisession, workers = 3)
progressr::handlers("progress")
results <- progressr::with_progress(run_nl_all(nl))

# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Write output to outpath of experiment within nl
write_simoutput(nl)

# Do further analysis:
analyze_nl(nl)