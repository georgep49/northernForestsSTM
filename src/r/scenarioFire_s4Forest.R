library(tidyverse)
library(vroom)
library(patchwork)
library(svglite)

source("src/r/wranglers.r")

# define types to deal with zero fire files....
ct <- c("idccliddddddddddddddddddddddddddi")
f <- list.files("src/data/s4/s4Forest/fireRecordsS4f", pattern = "fire_record", full.names = TRUE) |>
  str_sort(numeric = TRUE)

# this works so we have data.frames
X <- lapply(f, FUN = vroom, col_types = ct)

file_table <- tibble(file_name = f) |>
  mutate(ensemble = row_number())

fire_allreps <- bind_rows(X, .id = "ensemble") |>
  mutate(ensemble = as.numeric(ensemble)) |>
  left_join(file_table, by = "ensemble")

write_csv(fire_allreps, 
  file = "src/data/s4/s4Forest/fireDispersalLHC_allfire_records.csv")

#########
library(tidyverse)

class_names <- c("prop_gr", "prop_dSh", "prop_mSh", "prop_kshK", "prop_kshNoK", "prop_yfK", " prop_yfNok", "prop_old" , "prop_kshP", "prop_yfP", "prop_oldP")
class_names_topo <- c(paste0(class_names, "_gly"), paste0(class_names, "_slp"), paste0(class_names, "_rdg"))

names_lu <- read_csv("src/r/stateNames.csv")

state_fire_s4f_raw <- read_csv("src/data/s4/s4Forest/scenario4_forest_lhs.csv") |>
    janitor::clean_names() |>
    arrange(siminputrow, step)

fireDispersalLHC_allreps <- read_csv(file = "src/data/s4/s4Forest/fireDispersalLHC_allfire_records.csv")

fireLHC_size <- fireDispersalLHC_allreps |>
  group_by(ensemble) |>
  summarise(n = n(),
            total_size = sum(fire_size, na.rm = TRUE),
            mean_size = mean(fire_size, na.rm = TRUE),
            median_size = median(fire_size, na.rm = TRUE),
            max_size = max(fire_size, na.rm = TRUE)) |>
  ungroup()

# 1 - gully, 2 - slope, 3 - ridge
state_fire_s4f <- parse_to_state(state_fire_s4f_raw)

save.image("src/data/s4/s4Forest/s4ForestAllData.RData")
###

traps <- state_fire_s4f |>
  select(siminputrow, step, invasion, fire_frequency, seed_pred, flamm_start, extrinsic_sd, enso_freq_wgt,  farm_edge, farm_edge, terrain_type, starts_with("prop_"))

# get most prevalent type at end of run
# flat

traps_pal <- c("prop_dSh" = "#e7298a", "prop_mSh" = "#d95f02",
               "prop_kshNoK" = "#a6761d", "prop_kshK" = "#7570b3",
               "prop_yfK" = "#66a61e",  "prop_old" = "#1b9e77")

traps_flat <- traps |>
  filter(terrain_type == "flat")
dom <- apply(traps_flat[,11:21], 1, function(x) which(x == max(x)))
traps_flat$dom_state <- names(traps_flat[,11:21])[dom]
traps_flat$dom_abund <- apply(traps_flat[,11:21], 1, max) / (256 ^ 2)

traps_flat_gg <- ggplot(data = traps_flat, aes(x  = fire_frequency, y = seed_pred) ) +
  geom_point(aes(size = dom_abund, col = dom_state), alpha = 0.6) +
  scale_colour_manual(values = traps_pal) +
  scale_size_continuous(limits = c(0, 1), breaks = seq(0., 1.0, 0.2)) +
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
traps_ridge$dom_state <- names(traps_ridge[,11:21])[dom]
traps_ridge$dom_abund <- apply(traps_ridge[,11:21], 1, max) / (256 ^ 2)

traps_ridge_gg <- ggplot(data = traps_ridge, aes(x  = fire_frequency, y = seed_pred) ) +
  geom_point(aes(size = dom_abund, col = dom_state), alpha = 0.6) +
  scale_colour_manual(values = traps_pal) +   
  scale_size_continuous(limits = c(0, 1), breaks = seq(0., 1.0, 0.2)) +
  ggh4x::facet_nested_wrap(farm_edge ~ invasion, ncol = 1, nest_line =  TRUE, strip.position = "right") +   
  theme_bw() +
  theme(legend.position = "bottom", 
    strip.background = element_rect(fill = NA, color = NA), 
    ggh4x.facet.nestline = element_line(linetype = 3))

traps_gg_s4f <- traps_flat_gg | traps_ridge_gg +
    plot_annotation(tag_levels = "a") +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin())

library(svglite)
svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/figSM-fireDispersalLHCtrapsForest.svg", height = 13, width = 8, fix_text_size = FALSE)
traps_gg_s4f
dev.off()

save.image("src/data/s4/s4Forest/s4ForestAllData.RData")
