library(nlrx)

# This needs NetLogo 6.3.0 or earlier... to do with paths to extensions
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.3.0")
#modelpath <- file.path("H://files.fos/private/Research/NgaRakauTaketake/stm/src/nlogo/stmNorthernForests_63.nlogo")
modelpath <- file.path("e://Research/NgaRakauTaketake/stm/src/nlogo/stmNorthernForests_63.nlogo")
outpath <- file.path("E:/dump")

#######################
# SCENARIO 3 IN THE ECOGRAPHY PAPER
#######################

# Setup nl object
nl <- nl(nlversion = "6.3.0",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

# Attach experiment
nl@experiment <- experiment(expname = "test",
                            idrunnum = "nlrx-info",
                            outpath = outpath,
                            repetition = 1,
                            idsetup = "setup",
                            idgo = "go",
                            idfinal = "write-fire-record",
                            runtime = 10,
                            tickmetrics = "true",
                            metrics = c("abundances"),
                            variables = list("fire-frequency" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
                                             "flamm-start" = list(min = 0.0, max = 0.8, step = 0.01, qfun = "qunif"),
                                             "extrinsic-sd" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
                                             "enso-freq-wgt" = list(min = 0.9, max = 1.1, step = 0.01, qfun = "qunif"),
                                             "farm-edge?" = list(min = 0, max = 1, qfun = "qunif"),
                                             "invasion?" = list(min = 0, max = 1, qfun = "qunif"),
                                             "terrain-type" = list(min = 0, max = 1, qfun = "qunif"),
                                             "init-composition-file" = list(min = 0, max = 1, qfun = "qunif")),
                            constants = list("perc-seed" = 0.57,
                                             "fraction-seed-ldd" = 0.15,
                                             "seed-pred" = 0.0,
                                             "track-stalled?" = "false",
                                             "mean-ldd" = 5.0,
                                             "fire-slow" = 2,
                                             "fire-invasion" = 0.10,
                                             "base-invasion" = 0.05,
                                             "base-seed-prod-of" = 8,
                                             "base-seed-prod-yf" = 4,
                                             "crit-density-yng" = 12,
                                             "crit-density-old" = 14,
                                             "max-ticks" = 10,
                                             "burn-in-regen" = 10,
                                             "max-forest" = 1.0,
                                             "write-record?" = "true",
                                             "rust-global-inf" = 0.0,
                                             "phy-global-inf" = 0.0,
                                             "phy-local-inf" = 0.0,
                                             "phy-radius-inf" = 0.0,
                                             "sap-herbivory" = 0.0,
                                             "farm-edge-nodes" = 30,
                                             "mean-farm-depth" = 10,
                                             "farm-revegetate?" = "false"))


# Attach simdesign
nl@simdesign <- simdesign_lhs(nl = nl, samples = 2.5e3, nseeds = 1, precision = 3)

# apply transformation to boolean (netlogo needs booleans as strings)
# https://stackoverflow.com/questions/71067047/nlrx-package-include-boolean-parameter-as-variable
# note escape character around strings
nl@simdesign@siminput <- nl@simdesign@siminput %>% 
  dplyr::mutate(`farm-edge?` = dplyr::if_else(`farm-edge?` < 0.5, "false", "true")) |>
  dplyr::mutate(`invasion?` = dplyr::if_else(`invasion?` < 0.5, "false", "true")) |>
  dplyr::mutate(`terrain-type` = dplyr::if_else(`terrain-type` < 0.5, "\"flat\"", "\"ridge-gully\"")) |>
  dplyr::mutate(`forest-gully-prop` = dplyr::if_else(`terrain-type` == "\"flat\"", 0, 0.5))  |>
  dplyr::mutate(`init-composition-file` = dplyr::if_else(`init-composition-file` < 0.5,
                             "\"parameter_files/initial_shrub_composition.dat\"",
                             "\"parameter_files/initial_forest_composition.dat\""))

  
# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

results <- run_nl_all(nl)

# Run all simulations (loop over all siminputrows and simseeds)
library(future)
plan(multisession, workers = 14)
progressr::handlers("progress")
results <- progressr::with_progress(run_nl_all(nl))

# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Write output to outpath of experiment within nl
write_simoutput(nl)

# Do further analysis:
analyze_nl(nl)


################
# SCENARIO 4 IN THE ECOGRAPHY PAPER
################
library(nlrx)

# This needs NetLogo 6.3.0 or earlier... to do with paths to extensions
# Windows default NetLogo installation path (adjust to your needs!):
netlogopath <- file.path("C:/Program Files/NetLogo 6.3.0")
modelpath <- file.path("H://files.fos/private/Research/NgaRakauTaketake/stm/src/nlogo/stmNorthernForests_63.nlogo")
# modelpath <- file.path("D://Research/NgaRakauTaketake/stm/src/nlogo/stmNorthernForests_63.nlogo")
outpath <- file.path("C:/Temp")


# Setup nl object
nl <- nl(nlversion = "6.3.0",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

# Attach experiment
nl@experiment <- experiment(expname = "fireDispersal",
                            idrunnum = "nlrx-info",
                            outpath = outpath,
                            repetition = 1,
                            idsetup = "setup",
                            idgo = "go",
                            idfinal = "write-fire-record",
                            runtime = 300,
                            tickmetrics = "true",
                            metrics = c("abundances"),
                            variables = list("fire-frequency" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
                                             "flamm-start" = list(min = 0.0, max = 0.8, step = 0.01, qfun = "qunif"),
                                             "extrinsic-sd" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
                                             "enso-freq-wgt" = list(min = 0.9, max = 1.1, step = 0.01, qfun = "qunif"),
                                             "fraction-seed-ldd" = list(min = 0, max = 0.4, step = 0.02, qfun = "qunif"),
                                             "seed-pred" = list(min = 0, max = 0.6, step = 0.02, qfun = "qunif"),
                                             "mean-ldd" = list(min = 1, max = 10, step = 1, qfun = "qunif"),
                                             "sap-herbivory" = list(min = 0.0, max = 0.5, step = 0.02, qfun = "qunif"),
                                             "farm-edge?" = list(min = 0, max = 1, qfun = "qunif"),
                                             "terrain-type" = list(min = 0, max = 1, qfun = "qunif"),
                                             "init-composition-file" = list(min = 0, max = 1, qfun = "qunif")),
                            constants = list("perc-seed" = 0.57,
                                             "track-stalled?" = "false",
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
                                             "write-record?" = "true",
                                             "rust-global-inf" = 0.0,
                                             "phy-global-inf" = 0.0,
                                             "phy-local-inf" = 0.0,
                                             "phy-radius-inf" = 0.0,
                                             "farm-edge-nodes" = 30,
                                             "mean-farm-depth" = 10,
                                             "farm-revegetate?" = "false"))


# Attach simdesign
nl@simdesign <- simdesign_lhs(nl = nl, samples = 5e3, nseeds = 1, precision = 3)

# apply transformation to boolean (netlogo needs booleans as strings)
# https://stackoverflow.com/questions/71067047/nlrx-package-include-boolean-parameter-as-variable
nl@simdesign@siminput <- nl@simdesign@siminput |>
  dplyr::mutate(`farm-edge?` = dplyr::if_else(`farm-edge?` < 0.5, "false", "true")) |>
  dplyr::mutate(`terrain-type` = dplyr::if_else(`terrain-type` < 0.5, "\"flat\"", "\"ridge-gully\"")) |>
  dplyr::mutate(`forest-gully-prop` = dplyr::if_else(`terrain-type` == "\"flat\"", 0, 0.5))  |>
  dplyr::mutate(`init-composition-file` = dplyr::if_else(`init-composition-file` < 0.5,
                             "\"parameter_files/initial_shrub_composition.dat\"",
                             "\"parameter_files/initial_forest_composition.dat\""))

# Evaluate nl object:
eval_variables_constants(nl)
print(nl)

# Run all simulations (loop over all siminputrows and simseeds)
library(future)
plan(multisession, workers = 14)
progressr::handlers("progress")
results <- progressr::with_progress(run_nl_all(nl))

# Attach results to nl object:
setsim(nl, "simoutput") <- results

# Write output to outpath of experiment within nl
write_simoutput(nl)
