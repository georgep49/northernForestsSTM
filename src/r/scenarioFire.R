library(tidyverse)
library(vroom)

# define types to deal with zero fire files....
ct <- c("idccliddddddddddddddddddddddddddi")
f <- list.files("src/data/fireFarms/", pattern = "fire_record", full.names = TRUE) |>
  str_sort(numeric = TRUE)

# maybe this works so we have data.frames
X <- lapply(f, FUN = read_csv, col_types = ct)

file_table <- tibble(file_name = f) |>
  mutate(ensemble = row_number())

fire_allreps <- bind_rows(X, .id = "ensemble") |>
  mutate(ensemble = as.numeric(ensemble)) |>
  left_join(file_table, by = "ensemble")

library(poweRlaw)
m_pl <- displ$new(ff)
est <- estimate_xmin(m_pl)
m_pl$setXmin(est)
