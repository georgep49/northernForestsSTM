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
