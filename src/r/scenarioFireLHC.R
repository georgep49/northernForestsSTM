library(tidyverse)
library(vroom)

# define types to deal with zero fire files....
ct <- c("idccliddddddddddddddddddddddddddi")
f <- list.files("src/data/fireLHC/", pattern = "fire_record", full.names = TRUE) |>
  str_sort(numeric = TRUE)

# maybe this works so we have data.frames
X <- lapply(f, FUN = read_csv, col_types = ct)

file_table <- tibble(file_name = f) |>
  mutate(ensemble = row_number())

fire_allreps <- bind_rows(X, .id = "ensemble") |>
  mutate(ensemble = as.numeric(ensemble)) |>
  left_join(file_table, by = "ensemble")

write_csv(fire_allreps, file = "src/data/fireLHC/fireLHC_allfire_records.csv")


#########
library(tidyverse)

# variables = list("fire-frequency" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "flamm-start" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "extrinsic-sd" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "enso-freq-wgt" = list(min = 0.9, max = 1.1, step = 0.01, qfun = "qunif"),
#                  "farm-edge?" = list(min = 0, max = 1, qfun = "qunif"))

class_names <- c("prop_gr", "prop_dSh", "prop_mSh", "prop_kshK", "prop_kshNok", "prop_yfK", "prop_yfNok", "prop_old" , "prop_kshP", "prop_yfP", "prop_oldP")

names_lu <- read_csv("src/r/stateNames.csv")

state_fireLHC <- read_csv("src/data/fireLHC/fireClimate_lhc.csv") |>
    janitor::clean_names() |>
    arrange(siminputrow, step)

fireLHC_allreps <- read_csv(file = "src/data/fireLHC/fireLHC_allfire_records.csv")

fireLHC_size <- fireLHC_allreps |>
  group_by(ensemble) |>
  summarise(mean_size = mean(fire_size, na.rm = TRUE),
            median_size = median(fire_size, na.rm = TRUE),
            max_size = max(fire_size, na.rm = TRUE)) |>
  ungroup()

state_fireLHC <- state_fireLHC |>
    mutate(abundances = mgsub::mgsub(abundances, pattern = c("\\[", "\\]"), replacement = c("", ""))) |>
    separate_wider_delim(abundances, delim = " ", names = class_names) |>
    mutate(across(starts_with("prop_"), ~as.numeric(.x)))

state_fireLHC <- state_fireLHC |>
  left_join(fireLHC_size, by = c("siminputrow" = "ensemble")) |>
  rowwise() |>
  mutate(
    prop_ksh = sum(prop_kshK, prop_kshNok, prop_kshP),
    prop_yfor = sum(prop_yfK, prop_yfNok, prop_yfP),
    prop_ofor = sum(prop_old, prop_oldP)) |>
  ungroup()

summary(state_fireLHC$mean_size)


###

traps <- state_fireLHC |>
  select(siminputrow, step, invasion, fire_frequency, flamm_start, extrinsic_sd, enso_freq_wgt,  farm_edge, farm_edge, starts_with("prop_")) |>
  select(-(20:22)) |>
  filter(step == 300)


dom <- apply(traps[,9:19], 1, function(x) which(x == max(x)))
traps$dom_state <- names(traps[,9:19])[dom]
traps$dom_abund <- apply(traps[, 9:19], 1, max) / (256 ^ 2)

trap_gg <- ggplot(data = traps, aes(x  = fire_frequency, y = extrinsic_sd) ) +
  geom_point(aes(size = dom_abund, col = dom_state), alpha = 0.6) +
  facet_grid(farm_edge ~ invasion) +
  scale_colour_brewer(type = "qual", palette = "Dark2", direction = -1) +
  theme_bw()

library(svglite)
svglite(file = "figX-fireLHCtraps.svg", height = 8, width = 8, fix_text_size = FALSE)
trap_gg
dev.off()

####
library(tidyverse)
library(svglite)

fireHistory <- read_csv("src/data/fireLHC/fireLHC_allfire_records.csv") |>
  janitor::clean_names()

fireHistory <- mutate(fireHistory, fire_size_prop = fire_size / (256 ^ 2))

clim_fire <- ggplot(fireHistory %>% slice_sample(prop = 0.2)) +
      geom_point(aes(x = extrinsic, y = fire_size_prop, col = pre_prop_old_f)) +
      geom_vline(xintercept = 0) +
      facet_wrap(~start_farm) +
      scale_colour_distiller(palette = "Greens", direction = 1)

sf <- tibble(s = fireHistory$fire_size) |>
  arrange(-s) |>
  mutate(f = row_number())

sf_fire <- ggplot(data = sf) +
  geom_point(aes(x = f, y = s)) +
  scale_x_log10() +
  scale_y_log10()
 