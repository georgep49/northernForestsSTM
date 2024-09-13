
library(tidyverse)

class_names <- c("prop_gr", "prop_dSh", "prop_mSh", "prop_kshK", "prop_kshNok", "prop_yfK", "prop_yfNok", "prop_old" , "prop_kshP", "prop_yfP", "prop_oldP")

baseline <- read_csv("src/data/stmNorthernForests baseline-table.csv", skip = 6) |>
    janitor::clean_names() |>
    arrange(run_number, step)

baseline <- baseline |>
    mutate(abundances = mgsub::mgsub(abundances, pattern = c("\\[", "\\]"), replacement = c("", ""))) |>
    separate_wider_delim(abundances, delim = " ", names = class_names) |>
    mutate(across(starts_with("prop_"), ~as.numeric(.x)))


x <- baseline |>
    filter(step == 300) |>
    select(run_number, step, starts_with("prop_")) |>
    pivot_longer(-(1:2)) |>
    mutate(value = value / (256 ^ 2))

ggplot(x) +
    geom_violin(aes(x = name, y = value))



