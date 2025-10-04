library(tidyverse)
library(vroom)

# Start with the 'young' landscape runs
# define types to deal with zero fire files....
ct <- c("idccliddddddddddddddddddddddddddi")
f <- list.files("src/data/s2/fireRecordsScenario2", pattern = "fire_record", full.names = TRUE) |>
  str_sort(numeric = TRUE)

# works so we have data.frames
X <- lapply(f, FUN = vroom, col_types = ct)

# make sure we get the correct row numbers for netlogo join, so strip from the filename
file_table <- tibble(file_name = f) |>
  rowwise() |>
  mutate(ensemble = gsub(".csv", "", str_split(file_name, "_")[[1]][3])) |>
  mutate(ensemble = as.numeric(ensemble))

# again, check rows are lining up
fire_allreps <- bind_rows(X, .id = "rn") |>
  mutate(ensemble = as.numeric(file_table[rn,]$ensemble)) |>
  left_join(file_table, by = "ensemble")

vroom_write(fire_allreps, delim = ",", file = "src/data/s2/allfireRecords_s2.csv")

########
fire_records <- vroom(file = "src/data/s2/allfireRecords_s2.csv") |>
  janitor::clean_names()

fire_size <- fire_records |>
  group_by(ensemble) |>
  summarise(n = n(),
            total_size = sum(fire_size, na.rm = TRUE),
            mean_size = mean(fire_size, na.rm = TRUE),
            median_size = median(fire_size, na.rm = TRUE),
            max_size = max(fire_size, na.rm = TRUE)) |>
  ungroup()

####
class_names <- c("prop_gr", "prop_dSh", "prop_mSh", "prop_kshK", "prop_kshNok", "prop_yfK", "prop_yfNok", "prop_old" , "prop_kshP", "prop_yfP", "prop_oldP")
class_names_topo <- c(paste0(class_names, "_gly"), paste0(class_names, "_slp"), paste0(class_names, "_rdg"))

names_lu <- read_csv("src/r/stateNames.csv")

fire_farms <- vroom("src/data/s2/stmNorthernForests fire-invasion-forest-start-table.csv", skip = 6) |>
    janitor::clean_names() |>
    arrange(run_number, step)

# (slow)
fire_farms <- fire_farms |>
    mutate(abundances = str_remove_all(abundances, "\\[|\\]")) |>
    separate_wider_delim(abundances, delim = " ", names = class_names) |>
    mutate(text_data = str_remove_all(abundances_by_topo, "\\[|\\]")) |>  # remove brackets
    separate_wider_delim(cols = text_data, names = class_names_topo, delim = " ", too_few = "align_start") |>
    mutate(across(starts_with("prop_"), ~as.numeric(.x))) |>
  left_join(fire_size, by = c("run_number" = "ensemble")) |>
  rowwise() |>
  mutate(
    prop_ksh = sum(prop_kshK, prop_kshNok, prop_kshP),
    prop_yfor = sum(prop_yfK, prop_yfNok, prop_yfP),
    prop_ofor = sum(prop_old, prop_oldP))

state_fire_s2 <- filter(fire_farms, step == 300)

vroom_write(file = "src/data/s2/fireFarmsState_s2.csv", delim = ",", x = state_fire_s2)
save.image("src/data/s2/s2AllData.RData")

####################
lsp_type <- c(`parameter_files/initial_forest_composition.dat` = "Forest",
              `parameter_files/initial_shrub_composition.dat` = "Shrubland")

firesize_gg <- ggplot(data = state_fire_s2) +
  geom_jitter(aes(x = fire_frequency, y = mean_size / (256^2), col = farm_edge), 
                  alpha = 0.3, width = 0.0025) +
  geom_smooth(aes(x = fire_frequency, y = mean_size / (256^2), col = farm_edge))  +
  labs(x = "Fre frequency", y = "Mean fire size (prop. landscape)") +
  scale_colour_brewer(type = "qual", palette = "Dark2", name = "Farm edge?") +
  ggh4x::facet_nested_wrap(init_composition_file ~ terrain_type, 
        ncol = 1, nest_line =  TRUE, strip.position = "right",
        labeller = labeller(init_composition_file = lsp_type)) + 
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        ggh4x.facet.nestline = element_line(linetype = 3))

ofor_gg <- ggplot(data = state_fire_s2) +
  geom_jitter(aes(x = fire_frequency, y = prop_ofor / (256 ^ 2), col = farm_edge), alpha = 0.3, width = 0.0025) +
  geom_smooth(aes(x = fire_frequency, y = prop_ofor / (256 ^ 2), col = farm_edge))  +
  labs(x = "Fre frequency", y = "Proportional abundance old forest") +
  scale_colour_brewer(type = "qual", palette = "Dark2", name = "Farm edge?") +
  ggh4x::facet_nested_wrap(init_composition_file ~ terrain_type,
        ncol = 1, nest_line =  TRUE, strip.position = "left",
        labeller = labeller(init_composition_file = lsp_type)) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        ggh4x.facet.nestline = element_line(linetype = 3))


library(patchwork)
library(svglite)

freqFarm_gg <- ofor_gg + firesize_gg  +
  plot_layout(ncol = 2, guides = "collect", axes = "collect_x") +
  plot_annotation(tag_level = "a") &
  theme(legend.position = "bottom")

svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/fig3_freqFarm_plot.svg", 
  fix_text_size = FALSE, width = 8, height = 13)
freqFarm_gg
dev.off()

save.image("src/data/s2/s2AllData.RData")
################################

state_fire_s2_rg <- state_fire_s2 |>
  filter(terrain_type != "flat") |>
    mutate(tpi = str_remove_all(tpi_freq, "\\[|\\]")) |>
    separate_wider_delim(tpi, delim = " ", names = c("tpi_gully", "tpi_slope", "tpi_ridge")) |>
    mutate(across(starts_with("tpi_"), ~as.numeric(.x))) |>
    mutate(tpi_gully = tpi_gully / (256 ^ 2),
          tpi_slope = tpi_slope / (256 ^ 2),
          tpi_ridge = tpi_ridge / (256 ^ 2))

tpi_forest <- state_fire_s2_rg |>
  select(run_number, init_composition_file, terrain_type, farm_edge, fire_frequency, prop_old_gly, 
          prop_oldP_gly, prop_ofor, tpi_gully, tpi_slope, tpi_ridge) |>
  mutate(prop_ofor_gly = prop_old_gly + prop_oldP_gly) |>
  filter(prop_ofor > 0, init_composition_file == "parameter_files/initial_shrub_composition.dat") |>
  mutate(exp_gly_of = tpi_gully * prop_ofor) |>
  mutate(ratio_gly_of = prop_ofor_gly / exp_gly_of)
 
tpi_forest_gg <- ggplot(tpi_forest) + 
  geom_point(aes(x = prop_ofor / 256 ^ 2, y = ratio_gly_of, col = farm_edge), 
          alpha = 0.6, size = 2) +
  geom_hline(yintercept = 1, col = "red") +
  scale_colour_brewer(palette = "Dark2", name = "Farm edge?") +
  labs(x = "Proportion of old forest", y = "Expected:observed old forest") +
  theme_bw() +
  theme(legend.position = "bottom")

save.image("src/data/s2/s2AllData.RData")
