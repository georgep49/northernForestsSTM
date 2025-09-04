library(tidyverse)
library(vroom)

# define types to deal with zero fire files....
ct <- c("idccliddddddddddddddddddddddddddi")
f <- list.files("src/data/fireDispersalLHC/", pattern = "fire_record", full.names = TRUE) |>
  str_sort(numeric = TRUE)

# maybe this works so we have data.frames
X <- lapply(f, FUN = read_csv, col_types = ct)

file_table <- tibble(file_name = f) |>
  mutate(ensemble = row_number())

fire_allreps <- bind_rows(X, .id = "ensemble") |>
  mutate(ensemble = as.numeric(ensemble)) |>
  left_join(file_table, by = "ensemble")

write_csv(fire_allreps, file = "src/data/fireDispersalLHC/fireDispersalLHC_allfire_records.csv")


#########
library(tidyverse)

# variables = list("fire-frequency" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "flamm-start" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "extrinsic-sd" = list(min = 0.0, max = 0.2, step = 0.01, qfun = "qunif"),
#                  "enso-freq-wgt" = list(min = 0.9, max = 1.1, step = 0.01, qfun = "qunif"),
#                  "farm-edge?" = list(min = 0, max = 1, qfun = "qunif"))

class_names <- c("prop_gr", "prop_dSh", "prop_mSh", "prop_kshK", "prop_kshNok", "prop_yfK", "prop_yfNok", "prop_old" , "prop_kshP", "prop_yfP", "prop_oldP")

names_lu <- read_csv("src/r/stateNames.csv")

state_fireDispersalLHC <- read_csv("src/data/fireDispersalLHC/fireDispersal_lhc.csv") |>
    janitor::clean_names() |>
    arrange(siminputrow, step)

fireDispersalLHC_allreps <- read_csv(file = "src/data/fireDispersalLHC/fireDispersalLHC_allfire_records.csv")

fireDispersalLHC_size <- fireDispersalLHC_allreps |>
  group_by(ensemble) |>
  summarise(mean_size = mean(fire_size, na.rm = TRUE),
            median_size = median(fire_size, na.rm = TRUE),
            max_size = max(fire_size, na.rm = TRUE)) |>
  ungroup()

state_fireDispersalLHC <- state_fireDispersalLHC |>
    mutate(abundances = mgsub::mgsub(abundances, pattern = c("\\[", "\\]"), replacement = c("", ""))) |>
    separate_wider_delim(abundances, delim = " ", names = class_names) |>
    mutate(across(starts_with("prop_"), ~as.numeric(.x)))

state_fireDispersalLHC <- state_fireDispersalLHC |>
  left_join(fireDispersalLHC_size, by = c("siminputrow" = "ensemble")) |>
  rowwise() |>
  mutate(
    prop_ksh = sum(prop_kshK, prop_kshNok, prop_kshP),
    prop_yfor = sum(prop_yfK, prop_yfNok, prop_yfP),
    prop_ofor = sum(prop_old, prop_oldP)) |>
  ungroup()

summary(state_fireDispersalLHC$mean_size)


###

traps <- state_fireDispersalLHC |>
  select(siminputrow, step, invasion, fire_frequency, seed_pred, flamm_start, extrinsic_sd, enso_freq_wgt,  farm_edge, farm_edge, starts_with("prop_")) |>
  select(-(21:23)) |>
  filter(step == 300)


dom <- apply(traps[,10:20], 1, function(x) which(x == max(x)))
traps$dom_state <- names(traps[,10:20])[dom]
traps$dom_abund <- apply(traps[, 10:20], 1, max) / (256 ^ 2)

trap_gg <- ggplot(data = traps, aes(x  = fire_frequency, y = seed_pred) ) +
  geom_point(aes(size = dom_abund, col = dom_state), alpha = 0.6) +
  facet_grid(farm_edge ~ invasion) +
  scale_colour_manual(values = c("#E7298A","#7570B3", "#D95F02", "#1B9E77", "#66A61E")) +
  theme_bw()

library(svglite)
svglite(file = "figX-fireDispersalLHCtraps.svg", height = 8, width = 8, fix_text_size = FALSE)
trap_gg
dev.off()
