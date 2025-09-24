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
names_lu <- read_csv("src/r/stateNames.csv")


fire_farms <- vroom("src/data/s2/stmNorthernForests fire-invasion-forest-start-table.csv", skip = 6) |>
    janitor::clean_names() |>
    arrange(run_number, step)

# (slow)
fire_farms <- fire_farms |>
    mutate(abundances = mgsub::mgsub(abundances, pattern = c("\\[", "\\]"), replacement = c("", ""))) |>
    separate_wider_delim(abundances, delim = " ", names = class_names) |>
    mutate(across(starts_with("prop_"), ~as.numeric(.x))) |>
  left_join(fire_size, by = c("run_number" = "ensemble")) |>
  rowwise() |>
  mutate(
    prop_ksh = sum(prop_kshK, prop_kshNok, prop_kshP),
    prop_yfor = sum(prop_yfK, prop_yfNok, prop_yfP),
    prop_ofor = sum(prop_old, prop_oldP))

fire_farms_state <- filter(fire_farms, step == 300)

vroom_write(file = "src/data/s2/fireFarmsState_s2.csv", delim = ",", x = fire_farms_state)

save.image("src/data/s2/s2AllData.RData")

####################
lsp_type <- c(`parameter_files/initial_forest_composition.dat` = "Forest",
              `parameter_files/initial_shrub_composition.dat` = "Shrubland")

firesize_gg <- ggplot(data = fire_farms_state) +
  geom_jitter(aes(x = fire_frequency, y = mean_size / (256^2), col = farm_edge), 
                  alpha = 0.3, width = 0.0025) +
  geom_smooth(aes(x = fire_frequency, y = mean_size / (256^2), col = farm_edge))  +
  labs(x = "Fre frequency", y = "Mean fire size (prop. landscape)") +
  scale_colour_brewer(type = "qual", palette = "Dark2", name = "Farm edge?") +
  ggh4x::facet_nested_wrap(init_composition_file ~ terrain_type, 
        ncol = 1, nest_line =  TRUE, strip.position = "left",
        labeller = labeller(init_composition_file = lsp_type)) + 
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        ggh4x.facet.nestline = element_line(linetype = 3))

ofor_gg <- ggplot(data = fire_farms_state) +
  geom_jitter(aes(x = fire_frequency, y = prop_ofor / (256 ^ 2), col = farm_edge), alpha = 0.3, width = 0.0025) +
  geom_smooth(aes(x = fire_frequency, y = prop_ofor / (256 ^ 2), col = farm_edge))  +
  labs(x = "Fre frequency", y = "Proportional abundance old forest") +
  scale_colour_brewer(type = "qual", palette = "Dark2", name = "Farm edge?") +
  ggh4x::facet_nested_wrap(init_composition_file ~ terrain_type, 
        ncol = 1, nest_line =  TRUE, strip.position = "right",
        labeller = labeller(init_composition_file = lsp_type)) +   
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = NA, color = NA),
        ggh4x.facet.nestline = element_line(linetype = 3))

# degr_gg <- ggplot(data = fire_farms) +
#   geom_boxplot(aes(x = factor(fire_frequency), y = prop_dSh / (256^2), col = farm_edge)) +
#   labs(x = "Fre frequency", y = "Proportional abundance degraded") +
#   scale_colour_brewer(type = "qual", palette = "Dark2", name = "Farm edge?") +
#   theme_bw()

# manuka_gg <- ggplot(data = fire_farms) +
#   geom_boxplot(aes(x = factor(fire_frequency), y = prop_mSh / (256^2), col = farm_edge)) +
#   labs(x = "Fre frequency", y = "Proportional abundance manuka") +
#   scale_colour_brewer(type = "qual", palette = "Dark2", name = "Farm edge?") +
#   theme_bw()

library(patchwork)
library(svglite)

freqFarm_gg <- firesize_gg + ofor_gg +
  plot_layout(ncol = 2, guides = "collect", axes = "collect_x") +
  plot_annotation(tag_level = "a") &
  theme(legend.position = "bottom")

svglite(file = "../../Papers/Current/NSC/NRT/fire/figs/revised/fig3_freqFarm_plot.svg", fix_text_size = FALSE, width = 5.5, height = 8)
freqFarm_gg
dev.off()
####

save.image("src/data/s2/s2AllData.RData")
