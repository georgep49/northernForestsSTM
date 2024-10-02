#### Search for ranger model (allows search on ntrees unlike caret)
#### modified from https://bradleyboehmke.github.io/HOML/random-forest.html

library(tidyverse)
library(ranger)
library(pdp)
library(mlr)


ice_plot <- function(x, pred, cov_facet, cov_col)
{
  pred <- enquo(pred)
  cov_facet <- enquo(cov_facet)
  cov_col <- enquo(cov_col)

  ice_summ <- x %>%
    summarise(mn = mean(yhat), .by = c(!!pred, !!cov_facet) )

  if (quo_is_missing(cov_col))
  {
    g1 <- ggplot() +
      geom_line(data = x, aes(x = !!pred, y = yhat, group = yhat.id), col = "grey", alpha = 0.2) +
      geom_line(data = ice_summ, aes(x = !!pred, y = mn), col = "blue", linewidth = 1.5) +
      theme_bw()
  } else {
    g1 <- ggplot() +
      geom_line(data = x, aes(x = !!pred, y = yhat, group = yhat.id, col = !!cov_col), alpha = 0.2) +
      geom_line(data = ice_summ, aes(x = !!pred, y = mn), col = "blue", linewidth = 1.5) +
      scale_color_distiller(palette = "PiYG")
      theme_bw()
  }
  if (!quo_is_missing(cov_facet)) {g1 <- g1 + facet_wrap(vars(!!cov_facet))}

  list(ice_summ, g1)
}  


## Baseline in case we need it for reference
load("src/data/baseline/baseline.RData")

### Fire only SA
load("src/data/fireLHC/fireLHCRecords.RData")

fire_s2 <- state_fireLHC |>
  select(2:6, 33, 45, 46, 51) |>
  mutate(prop_ofor_ave = mean(prop_ofor / (256 ^2)),
    prop_ofor_min = min(prop_ofor) / (256 ^ 2),
   .by = siminputrow) |>
  filter(step == 300) |>
  mutate(prop_ofor = prop_ofor / (256 ^ 2),
      farm_edge = as.factor(farm_edge),
      mean_size = mean_size / (256 ^ 2))


## Old forest amount model
# Tune the hyperparameters using mlr and tuneRanger
fires2_of_task <- makeRegrTask(data = fire_s2[, c(1:5, 9)], target = "prop_ofor")
# estimateTimeTuneRanger(fires2_task)

# default measure is the variance reduction (~ impurity)
fires2_of_tune <- tuneRanger(fires2_of_task, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyperparameters
tuned_of_fires2 <- ranger(
  prop_ofor ~ .,
  data = fire_s2[, c(1:5, 9)],
  mtry = fires2_of_tune[[1]]$mtry,
  min.node.size = fires2_of_tune[[1]]$min.node.size,
  sample.fraction = fires2_of_tune[[1]]$sample.fraction,
  respect.unordered.factors = "order",
  importance = "impurity",
  seed = 314159)

## Fire size amount model
# Tune the hyperparameters using mlr and tuneRanger
fire_s2_compl <- fire_s2[complete.cases(fire_s2),]
fires2_area_task <- makeRegrTask(data = fire_s2_compl[, c(1:5, 8, 10)], target = "mean_size")
# estimateTimeTuneRanger(fires2_task)

# default measure is the variance reduction (~ impurity)
fires2_area_tune <- tuneRanger(fires2_of_task, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyperparameters
tuned_area_fires2 <- ranger(
  mean_size ~ .,
  data = fire_s2_compl[, c(1:5, 8, 10)],
  mtry = fires2_area_tune[[1]]$mtry,
  min.node.size = fires2_area_tune[[1]]$min.node.size,
  sample.fraction = fires2_area_tune[[1]]$sample.fraction,
  respect.unordered.factors = "order",
  importance = "impurity",
  seed = 314159)



################################
# Visualisation
####
library(iml)

# importance
predictor_of <- Predictor$new(model = tuned_of_fires2, data = fire_s2[,c(1:5)], y = fire_s2$prop_ofor)
imp_of <- FeatureImp$new(predictor_of, loss = "rmse")
imp_of_gg <- plot(imp_of) + 
  scale_x_continuous(limits = c(0, NA)) + 
  labs(x = "Feature importance", y = "Predictor") +
  theme_bw()

imp_of_gg

predictor_area <- Predictor$new(model = tuned_area_fires2, data = fire_s2_compl[,c(1:5, 10)], y = fire_s2_compl$mean_size)
imp_area <- FeatureImp$new(predictor_area, loss = "rmse")
imp_area_gg <- plot(imp_area) +
  scale_x_continuous(limits = c(0, NA)) + 
  labs(x = "Feature importance", y = "Predictor") +
  theme_bw()

imp_area_gg

# interactions
interact_all_of <- Interaction$new(predictor_of, grid.size = 15)
interact_of_gg <- plot(interact_all_of) +
  scale_x_continuous(limits = c(0, NA)) + 
  labs(x = "Feature importance", y = "Interaction strength") +
  theme_bw()

interact_of_gg

interact_all_area <- Interaction$new(predictor_area, grid.size = 15)
interact_area_gg <- plot(interact_all_area)  +
  scale_x_continuous(limits = c(0, NA)) + 
  labs(x = "Feature importance", y = "Interaction strength") +
  theme_bw()
interact_area_gg


#####
# ICE curves

p_firefreq <- pdp::partial(tuned_of_fires2, pred.var = "fire_frequency", ice = TRUE, plot.engine = "ggplot")
p_extrinsic <- pdp::partial(tuned_of_fires2, pred.var = "extrinsic_sd", ice = TRUE, plot.engine = "ggplot")

p_firefreq <- left_join(p_firefreq, fire_s2[,c(7, 3, 5)], by = c("yhat.id" = "siminputrow"))
p_extrinsic <- left_join(p_extrinsic, fire_s2[,c(7, 3, 5)], by = c("yhat.id" = "siminputrow"))

ice_firefreq <- ice_plot(p_firefreq, pred = fire_frequency, cov_facet = farm_edge, cov_col = extrinsic_sd)[[2]]
ice_extrinsic <- ice_plot(p_extrinsic, pred = extrinsic_sd, cov = farm_edge)[[2]]


library(patchwork)
library(svglite)
svglite(file = "scenario2_pdp.svg", width = 20, height = 16, fix_text_size = FALSE)
(p1 + ice_rgip) / (ice_rmp + ice_rsdp)
dev.off()


