library(tidyverse)
library(vroom)

# Start with the 'young' landscape runs
# define types to deal with zero fire files....
ct <- c("idccliddddddddddddddddddddddddddi")
f <- list.files("src/data/fireFarms/", pattern = "fire_record", full.names = TRUE) |>
  str_sort(numeric = TRUE)

# maybe this works so we have data.frames
X <- lapply(f, FUN = read_csv, col_types = ct)

file_table <- tibble(file_name = f) |>
  mutate(ensemble = row_number())

fireYng_allreps <- bind_rows(X, .id = "ensemble") |>
  mutate(ensemble = as.numeric(ensemble)) |>
  left_join(file_table, by = "ensemble")

write_csv(fire_allreps, file = "fireFarms_allfire_records.csv")
########

fireYng_records <- read_csv(file = "src/data/fireFarms/fireFarms_allfire_records.csv") |>
  janitor::clean_names()

fireYng_size <- fireYng_records |>
  group_by(ensemble) |>
  summarise(mean_size = mean(fire_size, na.rm = TRUE),
            median_size = median(fire_size, na.rm = TRUE),
            max_size = max(fire_size, na.rm = TRUE)) |>
  ungroup()


####
class_names <- c("prop_gr", "prop_dSh", "prop_mSh", "prop_kshK", "prop_kshNok", "prop_yfK", "prop_yfNok", "prop_old" , "prop_kshP", "prop_yfP", "prop_oldP")

names_lu <- read_csv("src/r/stateNames.csv")

fireYng_farms <- read_csv("src/data/fireFarms/stmNorthernForests fire-invasion-table.csv", skip = 6) |>
    janitor::clean_names() |>
    arrange(run_number, step)

fireYng_farms <- fireYng_farms |>
    mutate(abundances = mgsub::mgsub(abundances, pattern = c("\\[", "\\]"), replacement = c("", ""))) |>
    separate_wider_delim(abundances, delim = " ", names = class_names) |>
    mutate(across(starts_with("prop_"), ~as.numeric(.x)))

fireYng_farms <- fireYng_farms |>
  left_join(fireYng_size, by = c("run_number" = "ensemble")) |>
  rowwise() |>
  mutate(
    prop_ksh = sum(prop_kshK, prop_kshNok, prop_kshP),
    prop_yfor = sum(prop_yfK, prop_yfNok, prop_yfP),
    prop_ofor = sum(prop_old, prop_oldP))

########
# 'Old' forest runs now
ct <- c("idccliddddddddddddddddddddddddddi")
f <- list.files("src/data/fireFarmsForest/", pattern = "fire_record", full.names = TRUE) |>
  str_sort(numeric = TRUE)

# maybe this works so we have data.frames
X <- lapply(f, FUN = read_csv, col_types = ct)

file_table <- tibble(file_name = f) |>
  mutate(ensemble = row_number())

fireForest_allreps <- bind_rows(X, .id = "ensemble") |>
  mutate(ensemble = as.numeric(ensemble)) |>
  left_join(file_table, by = "ensemble")

write_csv(fireForest_allreps, file = "src/data/fireFarmsForest/fireFarmsForest_allfire_records.csv")
########

fireForest_records <- read_csv(file = "src/data/fireFarmsForest/fireFarmsForest_allfire_records.csv") |>
  janitor::clean_names()

fireForest_size <- fireForest_records |>
  group_by(ensemble) |>
  summarise(mean_size = mean(fire_size, na.rm = TRUE), 
            median_size = median(fire_size, na.rm = TRUE), 
            max_size = max(fire_size, na.rm = TRUE)) |>
  ungroup()

########

fireForest_farms <- read_csv("src/data/fireFarmsForest/stmNorthernForests fire-invasion-forest-start-table.csv", skip = 6) |>
    janitor::clean_names() |>
    arrange(run_number, step)

fireForest_farms <- fireForest_farms |>
    mutate(abundances = mgsub::mgsub(abundances, pattern = c("\\[", "\\]"), replacement = c("", ""))) |>
    separate_wider_delim(abundances, delim = " ", names = class_names) |>
    mutate(across(starts_with("prop_"), ~as.numeric(.x)))

fireForest_farms <- fireForest_farms |>
  left_join(fireForest_size, by = c("run_number" = "ensemble")) |>
  rowwise() |>
  mutate(
    prop_ksh = sum(prop_kshK, prop_kshNok, prop_kshP),
    prop_yfor = sum(prop_yfK, prop_yfNok, prop_yfP),
    prop_ofor = sum(prop_old, prop_oldP)) |>
  ungroup()

#######
# Merge the young and old starts....


fireForest_farms <- fireForest_farms |>
  mutate(run_number = run_number + 840,
         start_lsp = "old")

fireYng_farms <- fireYng_farms |>
         mutate(start_lsp = "yng")

fire_farms_state <- bind_rows(fireYng_farms, fireForest_farms %>% filter(step == 300))

write_csv(file = "src/data/fireFarms/fire_farms_state.csv", x = fire_farms_state)



####################
firesize_gg <- ggplot(data = fire_farms_state) +
  geom_jitter(aes(x = fire_frequency, y = mean_size / (256^2), col = farm_edge), alpha = 0.3, width = 0.005) +
  geom_smooth(aes(x = fire_frequency, y = mean_size / (256^2), col = farm_edge))  +
  labs(x = "Fre frequency", y = "Mean fire size (prop. landscape)") +
  scale_colour_brewer(type = "qual", palette = "Dark2", name = "Farm edge?") +
  facet_wrap(~start_lsp) +
  theme_bw() +
  theme(legend.position = "bottom")

ofor_gg <- ggplot(data = fire_farms_state) +
  geom_jitter(aes(x = fire_frequency, y = prop_ofor / (256 ^ 2), col = farm_edge), alpha = 0.3, width = 0.005) +
  geom_smooth(aes(x = fire_frequency, y = prop_ofor / (256 ^ 2), col = farm_edge))  +
  labs(x = "Fre frequency", y = "Proportional abundance old forest") +
  scale_colour_brewer(type = "qual", palette = "Dark2", name = "Farm edge?") +
  facet_wrap(~start_lsp) +
  theme_bw() +
  theme(legend.position = "bottom")

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
  plot_layout(nrow = 2, guides = "collect", axes = "collect_x") +
  plot_annotation(tag_level = "a") &
  theme(legend.position = "bottom")

svglite(file = "freqFarm_plot.svg", fix_text_size = FALSE, width = 5.5, height = 8)
freqFarm_gg
dev.off()
####


traps <- fire_farms_state |>
  select(fire_frequency, start_lsp, farm_edge, starts_with("prop_")) 

