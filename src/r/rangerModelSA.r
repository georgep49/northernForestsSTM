#### Search for ranger model (allows search on ntrees unlike caret)
#### modified from https://bradleyboehmke.github.io/HOML/random-forest.html

library(tidyverse)
library(ranger)
library(pdp)
library(mlr)
library(tuneRanger)
library(vivid)    # https://github.com/AlanInglis/vivid

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
      mean_size = ifelse(is.na(mean_size), 0, mean_size / (256 ^ 2)),
      gt_min_baseline = prop_ofor > 0.4504)  # 0.4504 is min ofor under baseline


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
fires2_area_task <- makeRegrTask(data = fire_s2[, c(1:5, 8, 10)], target = "mean_size")
# estimateTimeTuneRanger(fires2_task)

# default measure is the variance reduction (~ impurity)
fires2_area_tune <- tuneRanger(fires2_of_task, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyperparameters
tuned_area_fires2 <- ranger(
  mean_size ~ .,
  data = fire_s2[, c(1:5, 8, 10)],
  mtry = fires2_area_tune[[1]]$mtry,
  min.node.size = fires2_area_tune[[1]]$min.node.size,
  sample.fraction = fires2_area_tune[[1]]$sample.fraction,
  respect.unordered.factors = "order",
  importance = "impurity",
  seed = 314159)

## Better than min baseline amount model
# Tune the hyperparameters using mlr and tuneRanger
fires2_thresh_task <- makeClassifTask(data = fire_s2[, c(1:5, 12)], target = "gt_min_baseline")
#estimateTimeTuneRanger(fires2_thresh_task)

# default measure is the variance reduction (~ impurity)
fires2_thresh_tune <- tuneRanger(fires2_thresh_task, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyperparameters
tuned_thresh_fires2 <- ranger(
  gt_min_baseline ~ .,
  data = fire_s2[, c(1:5, 12)],
  classification = TRUE,
  mtry = fires2_thresh_tune[[1]]$mtry,
  min.node.size = fires2_thresh_tune[[1]]$min.node.size,
  sample.fraction = fires2_thresh_tune[[1]]$sample.fraction,
  respect.unordered.factors = "order",
  importance = "impurity",
  seed = 314159)


################################
# Visualisation
####
library(iml)
library(vivid)  # vivid from 

### mean fire size
vivi_area <- vivi(data = fire_s2[,c(1:5, 8, 10)],
                fit = tuned_area_fires2,
                response = "mean_size", 
                vars = c("prop_ofor_ave", "fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt") )

lbl_area <- data.frame(x = 1:6, y = 6:1, lbl = round(diag(vivi_area), 3))


viviHeat_area <- viviHeatmap(vivi_area) +
  geom_text(data = lbl_area, aes(x = x, y = y, label = lbl), size = 8) +
  theme(legend.position = "bottom")


pdp_area <- pdpVars(data = fire_s2[,c(1:5, 8, 10)],
         fit =  tuned_area_fires2, 
         response = "mean_size",      
         vars = c("prop_ofor_ave", "fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt"),
         nIce = 100)

pdp_area_trim <- pdp_area[[1]] + pdp_area[[2]] + pdp_area[[3]] + pdp_area[[4]] + pdp_area[[5]] + pdp_area[[6]]
  plot_layout(ncol = 3)

vivid_area <- viviHeat_area + pdp_area_trim +
  plot_layout(ncol = 2, widths = c(3,7))

####
### proportion of old forest after 300 y
vivi_of <- vivi(data = fire_s2[,c(1:5, 9)],
                fit = tuned_of_fires2,
                response = "prop_ofor", 
                vars = c("fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt") )

lbl_of <- data.frame(x = 1:5, y = 5:1, lbl = round(diag(vivi_of), 3))


viviHeat_of <- viviHeatmap(viviRf) +
  geom_text(data = lbl_of, aes(x = x, y = y, label = lbl), size = 8) +
  theme(legend.position = "bottom")


pdp_of <- pdpVars(data = fire_s2[,c(1:5, 9)],
         fit =  tuned_of_fires2, 
         response = "prop_ofor",      
         vars = c("fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt"),
         nIce = 100)

pdp_of_trim <- pdp_of[[1]] + pdp_of[[2]] + pdp_of[[3]] + pdp_of[[4]] + pdp_of[[5]] +
  plot_layout(ncol = 3)


vivid_of <- viviHeat_of + pdp_of_trim +
  plot_layout(ncol = 2, widths = c(3,7))

####
### threshold of min forest (TRUE class)
vivi_thresh <- vivi(data = fire_s2[,c(1:5,12)],
                fit = tuned_thresh_fires2,
                response = "gt_min_baseline",
                class = TRUE,
                vars = c("fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt") )

lbl_of <- data.frame(x = 1:5, y = 5:1, lbl = round(diag(vivi_thresh), 3))


viviHeat_thresh <- viviHeatmap(vivi_thresh) +
  geom_text(data = lbl_of, aes(x = x, y = y, label = lbl), size = 8) +
  theme(legend.position = "bottom")

####

library(svglite)
svglite(file = "Rplot%03d.svg", fix_text_size = FALSE, width = 12, height = 8)
vivid_area
vivid_of
viviHeat_thresh
dev.off()

######################################################################################################
# # importance
# predictor_of <- Predictor$new(model = tuned_of_fires2, data = fire_s2[,c(1:5)], y = fire_s2$prop_ofor)
# imp_of <- FeatureImp$new(predictor_of, loss = "rmse")
# imp_of_gg <- plot(imp_of) + 
#   scale_x_continuous(limits = c(0, NA)) + 
#   labs(x = "Feature importance", y = "Predictor") +
#   theme_bw()

# imp_of_gg

# ##
# predictor_area <- Predictor$new(model = tuned_area_fires2, data = fire_s2_compl[,c(1:5, 10)], y = fire_s2_compl$mean_size)
# imp_area <- FeatureImp$new(predictor_area, loss = "rmse")
# imp_area_gg <- plot(imp_area) +
#   scale_x_continuous(limits = c(0, NA)) + 
#   labs(x = "Feature importance", y = "Predictor") +
#   theme_bw()

# imp_area_gg

# ##
# predictor_thresh <- Predictor$new(model = tuned_thresh_fires2, data = fire_s2[,c(1:5)], y = fire_s2$gt_min_baseline)
# imp_thresh <- FeatureImp$new(predictor_thresh, loss = "ce")
# imp_thresh_gg <- plot(imp_thresh) +
#   scale_x_continuous(limits = c(0, NA)) + 
#   labs(x = "Feature importance", y = "Predictor") +
#   theme_bw()

# imp_thresh_gg



# # interactions
# interact_all_of <- Interaction$new(predictor_of, grid.size = 15)
# interact_of_gg <- plot(interact_all_of) +
#   scale_x_continuous(limits = c(0, NA)) + 
#   labs(x = "Feature importance", y = "Interaction strength") +
#   theme_bw()

# interact_of_gg

# ##
# interact_all_area <- Interaction$new(predictor_area, grid.size = 15)
# interact_area_gg <- plot(interact_all_area)  +
#   scale_x_continuous(limits = c(0, NA)) + 
#   labs(x = "Feature importance", y = "Interaction strength") +
#   theme_bw()
# interact_area_gg

# ##
# interact_all_thresh <- Interaction$new(predictor_thresh, grid.size = 15)
# interact_thresh_gg <- plot(interact_all_thresh)  +
#   scale_x_continuous(limits = c(0, NA)) + 
#   labs(x = "Feature importance", y = "Interaction strength") +
#   theme_bw()
# interact_thresh_gg


# library(patchwork)
# (imp_area_gg + imp_of_gg + imp_thresh_gg + interact_area_gg + interact_of_gg + interact_thresh_gg) +
#   plot_layout(nrow = 2) +
#   plot_annotation(tag_level = "a")



# #####
# # ICE curves

# # mean fire size

# p_ofor_ave_size <- pdp::partial(tuned_area_fires2, pred.var = "prop_ofor_ave", ice = TRUE, plot.engine = "ggplot")
# p_firefreq_size <- pdp::partial(tuned_area_fires2, pred.var = "fire_frequency", ice = TRUE, plot.engine = "ggplot")

# p_ofor_ave_size <- left_join(p_ofor_ave_size, fire_s2_compl[,c(7, 3, 5)], by = c("yhat.id" = "siminputrow"))
# p_firefreq_size <- left_join(p_firefreq_size, fire_s2_compl[,c(7, 3, 5)], by = c("yhat.id" = "siminputrow"))

# ice_ofor_ave_size <- ice_plot(p_ofor_ave_size, pred = prop_ofor_ave, cov_facet = farm_edge, cov_col = extrinsic_sd)[[2]]
# ice_firefreq_size <- ice_plot(p_firefreq_size, pred = fire_frequency, cov_facet = farm_edge, cov_col = extrinsic_sd)[[2]]

# x <- ice_ofor_ave_size + ice_firefreq_size
# x

# # ofor
# p_firefreq_of <- pdp::partial(tuned_of_fires2, pred.var = "fire_frequency", ice = TRUE, plot.engine = "ggplot")
# p_extrinsic_of <- pdp::partial(tuned_of_fires2, pred.var = "extrinsic_sd", ice = TRUE, plot.engine = "ggplot")

# p_firefreq_of <- left_join(p_firefreq_of, fire_s2[,c(7, 3, 5)], by = c("yhat.id" = "siminputrow"))
# p_extrinsic_of <- left_join(p_extrinsic_of, fire_s2[,c(7, 3, 5)], by = c("yhat.id" = "siminputrow"))

# ice_firefreq_of <- ice_plot(p_firefreq_of, pred = fire_frequency, cov_facet = farm_edge, cov_col = extrinsic_sd)[[2]]
# ice_extrinsic_of <- ice_plot(p_extrinsic_of, pred = extrinsic_sd, cov = farm_edge)[[2]]



# # min threshold
# p_firefreq_th <- pdp::partial(tuned_thresh_fires2, pred.var = "fire_frequency", ice = TRUE, plot.engine = "ggplot")
# p_extrinsic_th <- pdp::partial(tuned_thresh_fires2, pred.var = "extrinsic_sd", ice = TRUE, plot.engine = "ggplot")

# p_firefreq <- left_join(p_firefreq_thresh, fire_s2[,c(7, 3, 5)], by = c("yhat.id" = "siminputrow"))
# p_extrinsic <- left_join(p_extrinsic_thresh, fire_s2[,c(7, 3, 5)], by = c("yhat.id" = "siminputrow"))

# ice_firefreq <- ice_plot(p_firefreq_thresh, pred = fire_frequency, cov_facet = farm_edge, cov_col = extrinsic_sd)[[2]]
# ice_extrinsic <- ice_plot(p_extrinsic_thresh, pred = extrinsic_sd, cov = farm_edge)[[2]]