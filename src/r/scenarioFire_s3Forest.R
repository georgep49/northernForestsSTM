library(tidyverse)
library(vroom)
library(patchwork)
library(svglite)

source("src/r/wranglers.r")

# define types to deal with zero fire files....
ct <- c("idccliddddddddddddddddddddddddddi")
f <- list.files("src/data/s3/s3Forest/fireRecordsS3f", pattern = "fire_record", full.names = TRUE) |>
  str_sort(numeric = TRUE)

# maybe this works so we have data.frames
X <- lapply(f, FUN = vroom, col_types = ct)

file_table <- tibble(file_name = f) |>
  mutate(ensemble = row_number())

fire_allreps <- bind_rows(X, .id = "ensemble") |>
  mutate(ensemble = as.numeric(ensemble)) |>
  left_join(file_table, by = "ensemble")

write_csv(fire_allreps, file = "src/data/s3/s3Forest/fireLHC_allfire_records.csv")

#########
library(tidyverse)

# seems an issue with nlrx where it only records the abundances etc. as per the final tick??
# fine here as that is all we're analysing

# variables = list("fire-frequency" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "flamm-start" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "extrinsic-sd" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "enso-freq-wgt" = list(min = 0.9, max = 1.1, step = 0.01, qfun = "qunif"),
#                  "farm-edge?" = list(min = 0, max = 1, qfun = "qunif"))

class_names <- c("prop_gr", "prop_dSh", "prop_mSh", "prop_kshK", "prop_kshNok", 
  "prop_yfK", "prop_yfNok", "prop_old" , "prop_kshP", "prop_yfP", "prop_oldP")
class_names_topo <- c(paste0(class_names, "_gly"), paste0(class_names, "_slp"),
  paste0(class_names, "_rdg"))

names_lu <- read_csv("src/r/stateNames.csv")

state_fire_s3f_raw <- read_csv("src/data/s3/s3Forest/scenario3_forest_lhs.csv") |>
    janitor::clean_names() |>
    arrange(siminputrow, step)

fireLHC_allreps <- read_csv(file = "src/data/s3/s3Forest/fireLHC_allfire_records.csv")

fireLHC_size <- fireLHC_allreps |>
  group_by(ensemble) |>
  summarise(n = n(),
            total_size = sum(fire_size, na.rm = TRUE),
            mean_size = mean(fire_size, na.rm = TRUE),
            median_size = median(fire_size, na.rm = TRUE),
            max_size = max(fire_size, na.rm = TRUE)) |>
  ungroup()

# 1 - gully, 2 - slope, 3 - ridge
state_fire_s3f <- parse_to_state(state_fire_s3f_raw)

save.image("src/data/s3/s3Forest/s3ForestAllData.RData")

###
traps <- state_fire_s3f |>
  select(siminputrow, step, invasion, fire_frequency, flamm_start, 
          extrinsic_sd, enso_freq_wgt, farm_edge, farm_edge, terrain_type, 
          run_number, starts_with("prop_"))

# get most prevalent type at end of run
# flat

traps_pal <- c("prop_dSh" = "#e7298a", "prop_mSh" = "#d95f02", 
    "prop_kshK" = "#7570b3", "prop_yfK" = "#66a61e", "prop_old" = "#1b9e77")

traps_flat <- traps |> 
  filter(terrain_type == "flat")
dom <- apply(traps_flat[, 11:22], 1, function(x) which(x == max(x)))
traps_flat$dom_state <- names(traps_flat[, 11:22])[dom]
traps_flat$dom_abund <- apply(traps_flat[, 11:22], 1, max) / (256 ^ 2)

traps_flat_gg <- ggplot(data = traps_flat, aes(x  = fire_frequency, y = extrinsic_sd)) +
  geom_point(aes(size = dom_abund, col = dom_state), alpha = 0.6) +
  scale_colour_brewer(type = "qual", palette = "Dark2", direction = -1) +
  scale_size_continuous(limits = c(0, 1), breaks = seq(0.0, 1.0, 0.2)) +
  ggh4x::facet_nested_wrap(farm_edge ~ invasion, 
        ncol = 1, nest_line =  TRUE, strip.position = "left") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        ggh4x.facet.nestline = element_line(linetype = 3))
        
# ridge
traps_ridge <- traps |>
  filter(terrain_type != "flat")
dom <- apply(traps_ridge[,11:22], 1, function(x) which(x == max(x)))
traps_ridge$dom_state <- names(traps_ridge[, 11:22])[dom]
traps_ridge$dom_abund <- apply(traps_ridge[, 11:22], 1, max) / (256 ^ 2)

traps_ridge_gg <- ggplot(data = traps_ridge, aes(x  = fire_frequency, y = extrinsic_sd) ) +
  geom_point(aes(size = dom_abund, col = dom_state), alpha = 0.6) +
  scale_colour_brewer(type = "qual", palette = "Dark2", direction = -1) +
  scale_size_continuous(limits = c(0, 1), breaks = seq(0., 1.0, 0.2)) +
  ggh4x::facet_nested_wrap(farm_edge ~ invasion, 
        ncol = 1, nest_line =  TRUE, strip.position = "right") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        ggh4x.facet.nestline = element_line(linetype = 3))

traps_gg_s3f <- traps_flat_gg | traps_ridge_gg +
    plot_annotation(tag_levels = "a") +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin())


library(svglite)
svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/figSMX_fireLHCtraps_Forest.svg", height = 13, width = 8, fix_text_size = FALSE)
traps_gg_s3f
dev.off()

save.image("src/data/s3/s3Forest/s3ForestAllData.RData")
# load("src/data/s3/s3Forest/s3ForestAllData.RData")


####
library(tidyverse)
library(svglite)

fireHistory <- read_csv("src/data/s3Forest/fireLHC_allfire_records.csv") |>
  janitor::clean_names()

fireHistory <- mutate(fireHistory, fire_size_prop = fire_size / (256 ^ 2))

clim_fire <- ggplot(fireHistory %>% slice_sample(prop = 0.2)) +
      geom_point(aes(x = extrinsic, y = fire_size_prop, col = pre_prop_old_f)) +
      geom_vline(xintercept = 0) +
      facet_wrap(~start_farm) +
      scale_colour_distiller(palette = "Greens", direction = 1)

farmstart_fire <- ggplot(fireHistory %>% slice_sample(prop = 0.2)) +
  geom_histogram(aes(fire_size, after_stat(density))) +
  scale_x_log10() +
  facet_wrap(~start_farm)

save.image("src/data/s3/s3Forest/s3ForestAllData.RData")

########
# load("src/data/s3/s3Forest/s3ForestAllData.RData")

# gdata::keep(state_fire_s3f, sure = TRUE)

# tpi_frac <- vroom("src/data/baseline/stmNorthernForests tpi-fractions-table.csv", skip = 6) |>
#   janitor::clean_names()

# names(tpi_frac)[3:5] <- c("gully", "slope", "ridge")
# tpi_frac[3:5] <- tpi_frac[3:5] / (256 ^ 2)

# tpi_frac_summ <- tpi_frac |>
#   summarise(mean_gully = mean(gully), mean_slope = mean(slope), mean_ridge = mean(ridge))

# tpi_forest <- state_fire_s3f |>
#   select(run_number, terrain_type, farm_edge, fire_frequency, prop_old_gly, prop_oldP_gly, prop_ofor) |>
#   filter(terrain_type != "flat") |>
#   mutate(mean_gully = tpi_frac_summ$mean_gully) |>
#   mutate(prop_ofor_gly = prop_old_gly + prop_oldP_gly)

# x <- tpi_forest |>
#   filter(prop_ofor > 0) |>
#   mutate(exp_gly_of = mean_gully * prop_ofor) |>
#   mutate(ratio_gly_of = prop_ofor_gly / exp_gly_of)
 
#  ggplot(x) + 
#   geom_point(aes(x = prop_ofor, y = ratio_gly_of, col = farm_edge)) + 
#   geom_hline(yintercept = 1, col = "red") +
#   theme_bw()
