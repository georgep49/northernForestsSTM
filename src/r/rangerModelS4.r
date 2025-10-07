#### Search for ranger model (allows search on ntrees unlike caret)
#### modified from https://bradleyboehmke.github.io/HOML/random-forest.html

library(tidyverse)
library(ranger)
library(pdp)
library(mlr)
library(tuneRanger)
library(treeshap)
library(shapviz)

##################### 
# Shrublanb

### Fire only SA - merge forest and shrubland starts
load("src/data/s4/s4Shrub/s4ShrubAllData.RData")
load("src/data/s4/s4Forest/s4ForestAllData.RData")

state_fire_s4sf <- bind_rows(list(shrub = state_fire_s4s, forest = state_fire_s4f), .id = "start_lsp")

fire_s4sf <- state_fire_s4sf |>
  select(siminputrow, start_lsp, terrain_type, fire_frequency, flamm_start,  extrinsic_sd, fraction_seed_ldd, seed_pred, mean_ldd, sap_herbivory, enso_freq_wgt, farm_edge, invasion, step, prop_ofor, mean_size, terrain_type, start_lsp) |>
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

fire_s4sf <- split(fire_s4sf, f = list(fire_s4sf$start_lsp, fire_s4sf$terrain_type))

# titles for subsequent plots
ttl <- data.frame(sc = names(fire_s4sf)) |>
  mutate(sc = str_replace(sc, pattern = "\\.", replacement = "-")) |>
  mutate(sc = str_to_title(sc)) |>
  mutate(sc = str_replace(sc, pattern = "Ridge-", replacement = ""))


#############################
## Old forest amount model
ofor_params <- c("fire_frequency", "flamm_start",  "extrinsic_sd", "flamm_start", "seed_pred", "mean_ldd", "sap_herbivory", "enso_freq_wgt", "farm_edge", "invasion", "prop_ofor")

s4_of_rgr <- map(fire_s4sf, select, all_of(ofor_params))

# Tune the hyperparameters using mlr and tuneRanger (omit factors...)
# Use imap to iterate over the start_lsp x terrain combos
fires4_of_task <- imap(s4_of_rgr, ~ makeRegrTask(id = .y, data = .x, target = "prop_ofor"))
# estimateTimeTuneRanger(fires2_of_task)

# default measure is the variance reduction (~ impurity)
fires4_of_tune <- map(fires4_of_task, tuneRanger, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyper-parameters
rgr_of_hp <- map_dfr(fires4_of_tune, ~ {
    as.data.frame(as.data.frame(.x$recommended.pars))
  }, .id = "run")

# Run ranger row-wise
# Use pmap to iterate row-wise
# Convert parameter rows to list of lists (one per row)
rgr_of_hp <- transpose(rgr_of_hp)

# Run ranger row-wise with pmap
results_of <- pmap_dfr(
  list(params = rgr_of_hp, df = s4_of_rgr),
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
results_of <- results_of |>
  mutate(var_importance = map(model, ~ .x$variable.importance))

# results_of
# "forest.flat", "shrub.flat", "forest.ridge-gully", "shrub.ridge-gully"

vi_of <- map(results_of$model, pluck, "variable.importance") |>
  bind_rows() |>
  mutate(scenario = names(fires4_of_tune))  |>
  pivot_longer(cols = -scenario) |>
  group_by(scenario) |>
  mutate(value_sc = value / sum(value)) |>
  ungroup() |>
  separate_wider_delim(cols = scenario, delim = ".", 
        names = c("lsp", "terrain"), cols_remove = FALSE)

g_of <- ggplot(data = vi_of) +
  geom_col(aes(y = name, x = value_sc)) +
  ggh4x::facet_nested_wrap(lsp ~ terrain, 
        ncol = 1, nest_line =  TRUE, strip.position = "left") + 
  labs(x = "Scaled variable importance (impurity)", y = "Predictor (feature)") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        ggh4x.facet.nestline = element_line(linetype = 3)) 


# SHAP analysis and visualisation
s4_of_unify <- pmap(list(results_of$model, s4_of_rgr), unify)

# unified_rgr_of <- unify(results_of$model[[4]], s4_of_rgr[[4]])
s4_of_tX <- pmap(list(s4_of_unify, x = s4_of_rgr), treeshap)

s4_of_sv <- map(s4_of_tX, shapviz)
names(s4_of_sv) <- ttl$sc

# use the multi shapviz syntax 
# https://cran.csiro.au/web/packages/shapviz/vignettes/multiple_output.html

m_of <- mshapviz(s4_of_sv)
s4_of_vi <- sv_importance(m_of, show_numbers = TRUE, viridis_args = list())

s4_of_swarm <- sv_importance(m_of, kind = "beeswarm", bee_width = 0.25, show_numbers = TRUE) +
  plot_layout(ncol = 1)

s4_of_gg <- s4_of_vi / s4_of_swarm +
  plot_layout(heights = c(1, 4)) &
  theme_bw()

library(svglite)
svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/fig8-old_rangerS4.svg", height = 13, width = 8, fix_text_size = FALSE)
s4_of_gg
dev.off()

# patchwork::plot_layout(s4_of_svImpt, ncol = 1)

##################################
## Fire size amount model
# Tune the hyperparameters using mlr and tuneRanger

## Fire area model
area_params <- c("fire_frequency", "flamm_start",  "extrinsic_sd", "seed_pred", "mean_ldd", "sap_herbivory", "enso_freq_wgt", "farm_edge", "invasion", "mean_size")

s4_area_rgr <- map(fire_s4sf, select, all_of(area_params))

# Tune the hyperparameters using mlr and tuneRanger (omit factors...)
# USe imap to iterate over the start_lsp x terrain combos
fires4_area_task <- imap(s4_area_rgr, ~ makeRegrTask(id = .y, data = .x, target = "mean_size"))
# estimateTimeTuneRanger(fires4_area_task)

# default measure is the variance reduction (~ impurity)
fires4_area_tune <- map(fires4_area_task, tuneRanger, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyper-parameters
rgr_area_hp <- map_dfr(fires4_area_tune, ~ {
    as.data.frame(as.data.frame(.x$recommended.pars))
  }, .id = "run")

# Run ranger row-wise
# Use pmap to iterate row-wise
# Convert parameter rows to list of lists (one per row)
rgr_area_hp <- transpose(rgr_area_hp)

# Run ranger row-wise with pmap
results_area <- pmap_dfr(
  list(params = rgr_area_hp, df = s4_area_rgr),
  function(params, df) {
    model <- ranger(
      mean_size ~ .,
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
results_area <- results_area |>
  mutate(var_importance = map(model, ~ .x$variable.importance))


vi_area <- map(results_area$model, pluck, "variable.importance") |>
  bind_rows() |>
  mutate(scenario = names(fires4_area_tune))  |>
  pivot_longer(cols = -scenario) |>
  group_by(scenario) |>
  mutate(value_sc = value / sum(value)) |>
  ungroup() |>
  separate_wider_delim(cols = scenario, delim = ".", names = c("lsp", "terrain"), cols_remove = FALSE)

g_area <- ggplot(data = vi_area) +
  geom_col(aes(y = name, x = value_sc)) +
  ggh4x::facet_nested_wrap(lsp ~ terrain, 
        ncol = 1, nest_line =  TRUE, strip.position = "left") + 
  labs(x = "Scaled variable importance (impurity)", y = "Predictor (feature)") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        ggh4x.facet.nestline = element_line(linetype = 3))


# SHAP analysis and visualisation
s4_area_unify <- pmap(list(results_area$model, s4_area_rgr), unify)

# unified_rgr_of <- unify(results_of$model[[4]], s4_of_rgr[[4]])
s4_area_tX <- pmap(list(s4_area_unify, x = s4_area_rgr), treeshap)

s4_area_sv <- map(s4_area_tX, shapviz)
names(s4_area_sv) <- ttl$sc

# use the multi shapviz syntax 
# https://cran.csiro.au/web/packages/shapviz/vignettes/multiple_output.html

m_area <- mshapviz(s4_area_sv)
s4_area_vi <- sv_importance(m_area, show_numbers = TRUE, viridis_args = list())

s4_area_vi

s4_area_swarm <- sv_importance(m_area, kind = "beeswarm", bee_width = 0.25, show_numbers = TRUE) +
  plot_layout(ncol = 1)

s4_area_gg <- s4_area_vi / s4_area_swarm +
  plot_layout(heights = c(1, 4)) &
  theme_bw()

library(svglite)
svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/figSM-area_rangerS4.svg", 
  height = 13, width = 8, fix_text_size = FALSE)
s4_area_gg
dev.off()

save.image("src/data/s3s4Ranger/s4RangerModels.RData")

load("src/data/s3s4Ranger/s4RangerModels.RData")
