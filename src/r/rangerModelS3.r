#### Search for ranger model (allows search on ntrees unlike caret)
#### modified from https://bradleyboehmke.github.io/HOML/random-forest.html

library(tidyverse)
library(ranger)
library(pdp)
library(mlr)
library(tuneRanger)
library(treeshap)
library(shapviz)
#library(vivid)    # https://github.com/AlanInglis/vivid


##################### 
# Shrubland
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

### Fire only SA - merge forest and shrubland starts
load("src/data/s3/s3Shrub/s3ShrubAllData.RData")
load("src/data/s3/s3Forest/s3ForestAllData.RData")

state_fire_s3sf <- bind_rows(list(shrub = state_fire_s3s, forest = state_fire_s3f), 
  .id = "start_lsp")

fire_s3sf <- state_fire_s3sf |>
  select(siminputrow, start_lsp, terrain_type,fire_frequency, flamm_start, 
        extrinsic_sd, enso_freq_wgt, farm_edge, invasion, step, prop_ofor, delta_ofor,
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

# titles for subsequent plots
ttl <- data.frame(sc = names(fire_s3sf)) |>
  mutate(sc = str_replace(sc, pattern = "\\.", replacement = "-")) |>
  mutate(sc = str_to_title(sc)) |>
  mutate(sc = str_replace(sc, pattern = "Ridge-", replacement = ""))


#############################
## Old forest amount model
ofor_params <- c("fire_frequency", "flamm_start", "extrinsic_sd",  "enso_freq_wgt", "farm_edge", "invasion", "delta_ofor")

s3_of_rgr <- map(fire_s3sf, select, all_of(ofor_params))

# Tune the hyperparameters using mlr and tuneRanger (omit factors...)
# USe imap to iterate over the start_lsp x terrain combos
fires3_of_task <- imap(s3_of_rgr, ~ makeRegrTask(id = .y, data = .x, target = "delta_ofor"))
# estimateTimeTuneRanger(fires2_of_task)

# default measure is the variance reduction (~ impurity)
#  fires2_of_tune <- tuneRanger(fires2_of_task, num.trees = 1000, num.threads = 4, iters = 70)
fires3_of_tune <- map(fires3_of_task, tuneRanger, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyper-parameters
rgr_of_hp <- map_dfr(fires3_of_tune, ~ {
    as.data.frame(as.data.frame(.x$recommended.pars))
  }, .id = "run")

# Run ranger row-wise
# Use pmap to iterate row-wise
# Convert parameter rows to list of lists (one per row)
rgr_of_hp <- transpose(rgr_of_hp)

# Run ranger row-wise with pmap
results_of <- pmap_dfr(
  list(params = rgr_of_hp, df = s3_of_rgr),
  function(params, df) {
    model <- ranger(
      delta_ofor ~ .,
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
  mutate(scenario = names(fires3_of_tune))  |>
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
s3_of_unify <- pmap(list(results_of$model, s3_of_rgr), unify)

# unified_rgr_of <- unify(results_of$model[[4]], s3_of_rgr[[4]])
s3_of_tX <- pmap(list(s3_of_unify, x = s3_of_rgr), treeshap)

s3_of_sv <- map(s3_of_tX, shapviz)
names(s3_of_sv) <- ttl$sc

# use the multi shapviz syntax 
# https://cran.csiro.au/web/packages/shapviz/vignettes/multiple_output.html

m_of <- mshapviz(s3_of_sv)
s3_of_vi <- sv_importance(m_of, show_numbers = TRUE, viridis_args = list())

s3_of_swarm <- sv_importance(m_of, kind = "beeswarm", bee_width = 0.25, show_numbers = TRUE) +
  plot_layout(ncol = 1)

s3_of_gg <- s3_of_vi / s3_of_swarm +
  plot_layout(heights = c(1,4)) &
  theme_bw()


library(svglite)
svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/fig7-old_rangerS3_delta.svg", height = 13, width = 8, fix_text_size = FALSE)
s3_of_gg
dev.off()

##################################
## Fire size amount model
# Tune the hyperparameters using mlr and tuneRanger

## Fire area amount model
area_params <- c("fire_frequency", "flamm_start", "extrinsic_sd",  "enso_freq_wgt", "farm_edge", "invasion", "mean_size")

s3_area_rgr <- map(fire_s3sf, select, all_of(area_params))

# Tune the hyperparameters using mlr and tuneRanger (omit factors...)
# USe imap to iterate over the start_lsp x terrain combos
fires3_area_task <- imap(s3_area_rgr, ~ makeRegrTask(id = .y, data = .x, target = "mean_size"))
# estimateTimeTuneRanger(fires3_area_task)

# default measure is the variance reduction (~ impurity)
fires3_area_tune <- map(fires3_area_task, tuneRanger, num.trees = 1000, num.threads = 4, iters = 70)

# Build the old forest model with tuned hyper-parameters
rgr_area_hp <- map_dfr(fires3_area_tune, ~ {
    as.data.frame(as.data.frame(.x$recommended.pars))
  }, .id = "run")

# Run ranger row-wise
# Use pmap to iterate row-wise
# Convert parameter rows to list of lists (one per row)
rgr_area_hp <- transpose(rgr_area_hp)

# Run ranger row-wise with pmap
results_area <- pmap_dfr(
  list(params = rgr_area_hp, df = s3_area_rgr),
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
  mutate(scenario = names(fires3_area_tune))  |>
  pivot_longer(cols = -scenario) |>
  group_by(scenario) |>
  mutate(value_sc = value / sum(value)) |>
  ungroup() |>
  separate_wider_delim(cols = scenario, delim = ".", 
        names = c("lsp", "terrain"), cols_remove = FALSE)

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
s3_area_unify <- pmap(list(results_area$model, s3_area_rgr), unify)

# unified_rgr_of <- unify(results_of$model[[4]], s3_of_rgr[[4]])
s3_area_tX <- pmap(list(s3_area_unify, x = s3_area_rgr), treeshap)

s3_area_sv <- map(s3_area_tX, shapviz)
names(s3_area_sv) <- ttl$sc

# use the multi shapviz syntax 
# https://cran.csiro.au/web/packages/shapviz/vignettes/multiple_output.html

m_area <- mshapviz(s3_area_sv)
s3_area_vi <- sv_importance(m_area, show_numbers = TRUE, viridis_args = list())

s3_area_vi

s3_area_swarm <- sv_importance(m_area, kind = "beeswarm", bee_width = 0.25, show_numbers = TRUE) +
  plot_layout(ncol = 1)

s3_area_gg <- s3_area_vi / s3_area_swarm +
  plot_layout(heights = c(1,4)) &
  theme_bw()

library(svglite)
svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/figSM-area_rangerS3.svg", height = 13, width = 8, fix_text_size = FALSE)
s3_area_gg
dev.off()

save.image("src/data/s3s3Ranger/s3RangerModels.RData")
load("src/data/s3s3Ranger/s3RangerModels.RData")


#####
### Code to extract interactions
library(tidyverse)
library(shapviz)
load("src/data/s3s4Ranger/s3RangerModels.RData")

s3_of_shaps <- map(s3_of_sv, function(x) {colMeans(abs(x$S))} ) |> 
  bind_rows()

s3_of_shap_max <- colnames(s3_of_shaps)[apply(s3_of_shaps, 1, which.max)]

s3_of_ia_top <- pmap(list(s3_of_sv, s3_of_shap_max), potential_interactions) |> 
  bind_rows() |>
  mutate(expt = ttl$sc) |>
  separate_wider_delim(expt, delim = "-", names = c("lsp", "terrain"), cols_remove = FALSE) |>
  janitor::clean_names() |>
  pivot_longer(cols = -c(lsp, terrain, expt)) |>
  filter(!is.na(value))

var_tag <- data.frame(expt = ttl$sc, v_name = s3_of_shap_max) |>
  separate_wider_delim(expt, delim = "-", names = c("lsp", "terrain"), cols_remove = FALSE)


s3_of_ia_gg <- ggplot(s3_of_ia_top) +
  geom_col(aes(x = name, y = value, fill = value)) +
  facet_grid(lsp ~ terrain) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 0.8)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  geom_text(data = var_tag, aes(x = -Inf, y = Inf, label = paste("Main = ", v_name)), 
    hjust = -0.1, vjust = 1.5) +
  theme_bw() +
  theme()

#
s3_area_shaps <- map(s3_area_sv, function(x) {colMeans(abs(x$S))} ) |> 
  bind_rows()

s3_area_shap_max <- colnames(s3_area_shaps)[apply(s3_area_shaps, 1, which.max)]

s3_area_ia_top <- pmap(list(s3_area_sv, s3_area_shap_max), potential_interactions) |> 
  bind_rows() |>
  mutate(expt = ttl$sc) |>
  separate_wider_delim(expt, delim = "-", names = c("lsp", "terrain"), cols_remove = FALSE) |>
  janitor::clean_names() |>
  pivot_longer(cols = -c(lsp, terrain, expt)) |>
  filter(!is.na(value))

var_tag <- data.frame(expt = ttl$sc, v_name = s3_area_shap_max) |>
  separate_wider_delim(expt, delim = "-", names = c("lsp", "terrain"), cols_remove = FALSE)


s3_area_ia_gg <- ggplot(s3_area_ia_top) +
  geom_col(aes(x = name, y = value, fill = value)) +
  facet_grid(lsp ~ terrain) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 0.8)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  geom_text(data = var_tag, aes(x = -Inf, y = Inf, label = paste("Main = ", v_name)), 
    hjust = -0.1, vjust = 1.5) +
  theme_bw() +
  theme()

library(patchwork)
s3_inter_gg <- s3_of_ia_gg / s3_area_ia_gg +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", axis_titles = "collect")

library(svglite)
svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/SM/figSM-interactions_rangerS3.svg",
  height = 13, width = 8, fix_text_size = FALSE)
s3_inter_gg
dev.off()
