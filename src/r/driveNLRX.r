library(nlrx)

# This needs NetLogo 6.3.0 or earlier... to do with paths to extensions
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.3.0")
# modelpath <- file.path("H://files.fos/private/Research/NgaRakauTaketake/stm/src/nlogo/stmNorthernForests_63.nlogo")
modelpath <- file.path("D://Research/NgaRakauTaketake/stm/src/nlogo/stmNorthernForests_63.nlogo")
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
                            evalticks = "true",
                            metrics=c("abundances"),
                            variables = list("fire-frequency" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
                                             "flamm-start" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
                                             "extrinsic-sd" = list(min = -0.1, max = 0.1, step = 0.01, qfun = "qunif"),
                                             "enso-freq-wgt" = list(min = 0.9, max = 1.1, step = 0.01, qfun = "qunif")),
                            constants = list("perc-seed" = 0.57,
                                             "fraction-seed-ldd" = 0.15,
                                             "seed-pred" = 0.0,
                                             "track-stalled?" = FALSE,
                                             "mean-ldd" = 5.0,
                                             "base-invasion" = 0.05,
                                             "fire-slow" = 2,
                                             "fire-invasion" = 0.10,
                                             "base-seed-prod-of" = 8,
                                             "base-seed-prod-yf" = 4,
                                             "crit-density-yng" = 12,
                                             "crit-density-old" = 14,
                                             "max-ticks" = 300,
                                             "burn-in-regen" = 10,
                                             "max-forest" = 1.0,
                                             "write-record?" = TRUE,
                                             "extrinsic-sd" = 0.0,
                                             "rust-global-inf" = 0.0,
                                             "phy-global-inf" = 0.0,
                                             "phy-local-inf" = 0.0,
                                             "phy-radius-inf" = 1.5,
                                             "enso-matrix-file" = "",
                                             "sap-herbivory" = 0.0,
                                             "enso-freq-wgt" = 1.0,
                                             "farm-edge?" = FALSE,
                                             "farm-edge-nodes" = 30,
                                             "mean-farm-depth" = 10,
                                             "farm-revegetate?" = true,
                                             "nlrx-info" = ""))

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