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

##################### 
# Shrublanb
## Baseline in case we need it for reference
load("src/data/baseline/baseline.RData")

summ_old <- baseline |>
  group_by(start_lsp, terrain_type) |>
   summarise(
    min_of = min(prop_old + prop_oldP) / (256 ^ 2),
    mdn_of = median(prop_old + prop_oldP) / (256 ^ 2),
    mean_of = mean(prop_old + prop_oldP) / (256 ^ 2),
    max_of = max(prop_old + prop_oldP) / (256 ^ 2)
   ) |>
   ungroup()

### Fire only SA - shrubland start
load("src/data/s3Shrub/s3ShrubAllData.RData")
load("src/data/s3Forest/s3ForestAllData.RData")

state_fire_s3sf <- bind_rows(list(shrub = state_fire_s3s, forest = state_fire_s3f), 
  .id = "start_lsp")

fire_s3sf <- state_fire_s3sf |>
  select(siminputrow, start_lsp, terrain_type,fire_frequency, flamm_start, 
        extrinsic_sd, enso_freq_wgt, farm_edge, invasion, step, prop_ofor, 
        mean_size, terrain_type, start_lsp) |>
  mutate(prop_ofor_ave = mean(prop_ofor / (256 ^ 2))) |>
  group_by(siminputrow) |>
  mutate(siminputrow = paste0(siminputrow, str_sub(start_lsp, 1, 1)),
    prop_ofor_min = min(prop_ofor) / (256 ^ 2)) |>
  ungroup() |>
  filter(step == 300) |>
  mutate(prop_ofor = prop_ofor / (256 ^ 2),
      farm_edge = as.factor(farm_edge),
      invasion = as.factor(invasion),
      mean_size = ifelse(is.na(mean_size), 0, mean_size / (256 ^ 2)))

fire_s3sf <- split(fire_s3sf, f = list(fire_s3sf$start_lsp, fire_s3sf$terrain_type))

## Old forest amount model
ofor_params <- c("fire_frequency", "flamm_start", "extrinsic_sd",  "enso_freq_wgt", "farm_edge", "invasion", "prop_ofor")

s3_of_rgr <- map(fire_s3sf, select, all_of(ofor_params))

# Tune the hyperparameters using mlr and tuneRanger (omit factors...)
# USe imap to iterate over the start_lsp x terrain combos
fires3_of_task <- imap(s3_of_rgr, ~ makeRegrTask(id = .y, data = .x, target = "prop_ofor"))

# fires2_of_task <- makeRegrTask(data = fire_s3s$flat[, c(1:6, 8)], target = "prop_ofor")
# estimateTimeTuneRanger(fires2_of_task)

# default measure is the variance reduction (~ impurity)
#  fires2_of_tune <- tuneRanger(fires2_of_task, num.trees = 1000, num.threads = 4, iters = 70)
fires3_of_tune <- map(fires3_of_task, tuneRanger, num.trees = 1000, num.threads = 4, iters = 70)


# Build the old forest model with tuned hyper-parameters
rgr_hp <- map_dfr(fires3_of_tune, ~ {
    as.data.frame(as.data.frame(.x$recommended.pars))
  }, .id = "run")

# Run ranger row-wise
# Use pmap to iterate row-wise
# Convert parameter rows to list of lists (one per row)
rgr_hp <- transpose(rgr_hp)

# Run ranger row-wise with pmap
results <- pmap_dfr(
  list(params = rgr_hp, df = s3_of_rgr),
  function(params, df) {
    model <- ranger(
      prop_ofor ~ .,
      data = df,
      mtry = params$mtry,
      min.node.size = params$min.node.size,
      sample.fraction = params$sample.fraction,
      respect.unordered.factors = "order",
      importance = "impurity",
    )
    tibble(
      mtry = params$mtry,
      min.node.size = params$min.node.size,
      sample.fraction = params$sample.fraction,
      oob_error = model$prediction.error,
      model = list(model)
    )
  }
)

## add the variable importance
results <- results |>
  mutate(var_importance = map(model, ~ .x$variable.importance) )

######
# tuned_of_fires2 <- ranger(
#   prop_ofor ~ .,
#   data = fire_s3s$flat[, c(1:6, 8)],
#   mtry = fires2_of_tune[[1]]$mtry,
#   min.node.size = fires2_of_tune[[1]]$min.node.size,
#   sample.fraction = fires2_of_tune[[1]]$sample.fraction,
#   respect.unordered.factors = "order",
#   importance = "impurity",
#   seed = 314159)

## Fire size amount model
# Tune the hyperparameters using mlr and tuneRanger
fires2_area_task <- makeRegrTask(data = fire_s3[, c(1:6, 10, 11)], target = "mean_size")
# estimateTimeTuneRanger(fires2_task)

# default measure is the variance reduction (~ impurity)
fires2_area_tune <- tuneRanger(fires2_area_task, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyperparameters
tuned_area_fires2 <- ranger(
  mean_size ~ .,
  data = fire_s3[, c(1:6, 10, 11)],
  mtry = fires2_area_tune[[1]]$mtry,
  min.node.size = fires2_area_tune[[1]]$min.node.size,
  sample.fraction = fires2_area_tune[[1]]$sample.fraction,
  respect.unordered.factors = "order",
  importance = "impurity",
  seed = 314159)

## Better than min baseline amount model
# Tune the hyperparameters using mlr and tuneRanger
fires2_thresh_task <- makeClassifTask(data = fire_s3[, c(1:6, 13)], target = "gt_min_baseline")
#estimateTimeTuneRanger(fires2_thresh_task)

# default measure is the variance reduction (~ impurity)
fires2_thresh_tune <- tuneRanger(fires2_thresh_task, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyperparameters
tuned_thresh_fires2 <- ranger(
  gt_min_baseline ~ .,
  data = fire_s3[, c(1:6, 13)],
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
library(patchwork)
load("src/data/fireLHC/fireLHC.RData")

load("src/data/fireLHC/fireLHC_rangerModels.RData")

### mean fire size
vivi_area <- vivi(data = fire_s2[,c(1:6, 10, 11)],
                fit = tuned_area_fires2,
                response = "mean_size", 
                vars = c("prop_ofor_ave", "fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt", "invasion"))

lbl_area <- data.frame(x = 1:7, y = 7:1, lbl = round(diag(vivi_area), 3))


viviHeat_area <- viviHeatmap(vivi_area) +
  geom_text(data = lbl_area, aes(x = x, y = y, label = lbl), size = 8) +
  theme(legend.position = "bottom")


pdp_area <- pdpVars(data = fire_s2[,c(1:6, 10, 11)],
         fit =  tuned_area_fires2,
         response = "mean_size",
         vars = c("prop_ofor_ave", "fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt", "invasion"),
         nIce = 100)

pdp_area_top <- pdpVars(data = fire_s2[,c(1:6, 10, 11)],
                    fit =  tuned_area_fires2,
                    response = "mean_size",
                    vars = c("prop_ofor_ave", "extrinsic_sd", "farm_edge"),
                    nIce = 100)


# pdp_area_trim <- pdp_area[[1]] + pdp_area[[2]] + pdp_area[[3]] + pdp_area[[4]] + pdp_area[[5]] + pdp_area[[6]] + pdp_area[[7]] +
#  plot_layout(ncol = 4)

# using patchwork seems to drop the guide?
pdp_area_trim <- cowplot::plot_grid(plotlist = pdp_area_top)

vivid_area <- viviHeat_area + pdp_area_trim +
  plot_layout(ncol = 2, widths = c(3, 7))

####
### proportion of old forest after 300 y
vivi_of <- vivi(data = fire_s3[,c(1:6, 8)],
                fit = tuned_of_fires2,
                response = "prop_ofor", 
                vars = c("fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt", "invasion"))

lbl_of <- data.frame(x = 1:6, y = 6:1, lbl = round(diag(vivi_of), 3))


viviHeat_of <- viviHeatmap(vivi_of) +
  geom_text(data = lbl_of, aes(x = x, y = y, label = lbl), size = 8) +
  theme(legend.position = "bottom")


pdp_of <- pdpVars(data = fire_s3[,c(1:6, 8)],
         fit =  tuned_of_fires2, 
         response = "prop_ofor",      
         vars = c("fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt", "invasion"),
         nIce = 100)


pdp_of_top <- pdpVars(data = fire_s2[,c(1:6, 8)],
                  fit =  tuned_of_fires2, 
                  response = "prop_ofor",      
                  vars = c("fire_frequency", "extrinsic_sd", "farm_edge"),
                  nIce = 100)

# pdp_of_trim <- pdp_of[[1]] + pdp_of[[2]] + pdp_of[[3]] + pdp_of[[4]] + pdp_of[[5]] + pdp_of[[6]] +
#   plot_layout(ncol = 3)

pdp_of_trim <- cowplot::plot_grid(plotlist = pdp_of_top)

vivid_of <- viviHeat_of + pdp_of_trim +
  plot_layout(ncol = 2, widths = c(3,7))

####
### threshold of min forest (TRUE class)
vivi_thresh <- vivi(data = fire_s3[,c(1:6,13)],
                fit = tuned_thresh_fires2,
                response = "gt_min_baseline",
                class = TRUE,
                vars = c("fire_frequency", "extrinsic_sd", "farm_edge", "flamm_start", "enso_freq_wgt", "invasion") )

lbl_thresh <- data.frame(x = 1:6, y = 6:1, lbl = round(diag(vivi_thresh), 3))

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
