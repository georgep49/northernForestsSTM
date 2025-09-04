# Baseline plots

library(tidyverse)

class_names <- c("prop_gr", "prop_dSh", "prop_mSh", "prop_kshK", "prop_kshNok", "prop_yfK", "prop_yfNok", "prop_old" , "prop_kshP", "prop_yfP", "prop_oldP")

names_lu <- read_csv("src/r/stateNames.csv")

baseline_raw <- read_csv("src/data/baseline/stmNorthernForests baseline-shrub-ridge-table.csv", skip = 6) |>
    janitor::clean_names() |>
    arrange(run_number, step)

baseline <- baseline_raw |>
    mutate(abundances = mgsub::mgsub(abundances, pattern = c("\\[", "\\]"), replacement = c("", ""))) |>
    separate_wider_delim(abundances, delim = " ", names = class_names) |>
    mutate(across(starts_with("prop_"), ~as.numeric(.x))) |>
    mutate(start_lsp = if_else(str_detect(init_composition_file, "forest"), "forest", "shrubland"))

baseline_long <- baseline |>
    select(run_number, start_lsp, terrain_type, step, starts_with("prop_")) |>
    pivot_longer(cols = starts_with("prop_"), names_to = "state", values_to = "prop") |>
    mutate(prop = prop / (256^2)) |>
    left_join(names_lu, by = "state") |>
    filter(state != "prop_gr" & state != "prop_dSH")

time_states <- baseline_long |>
    group_by(step, state, start_lsp, terrain_type) |>
    summarise(as_tibble_row(quantile(prop, c(0.1, 0.5, 0.9)))) |>
    rename(prop10 = `10%`, median = `50%`, prop90 = `90%`) |>
    left_join(names_lu, by = "state")

final_state <- baseline_long |>
    filter(step == 300) |>
    mutate(state_label = forcats::fct_reorder(as.factor(state_label), prop, .desc = TRUE))


###
baseline_time_gg <- ggplot(data = time_states) +
    geom_line(aes(x = step, y = median, col = state_label)) +
    geom_ribbon(aes(x = step, ymin = prop10, ymax = prop90, fill = state_label), alpha = 0.3) + 
    ggrepel::geom_text_repel(data = time_states %>% filter(step == 300),
            aes(x = step, y = median, label = state_label),  na.rm = TRUE) +
    facet_grid(terrain_type ~ start_lsp) +
    theme_bw() +
    theme(legend.position = "bottom")


baseline_final_gg <- ggplot(data = final_state, aes(x = state_label, y = prop)) +
    geom_boxplot(aes(fill = state_group), outliers = FALSE) +
    geom_jitter(aes(fill = state_group, col = state_group), width = 0.1, alpha = 0.2) +
    labs(x = "State", y = "Final proportional abundance") +
    scale_fill_brewer(type = "qual", palette = "Dark2") +
    scale_colour_brewer(type = "qual", palette = "Dark2") +
    facet_grid(terrain_type ~ start_lsp) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1), 
            legend.position = "bottom")


library(patchwork)
library(svglite)

baseline_gg <- baseline_time_gg + baseline_final_gg +
  plot_annotation(tag_level = "a")

svglite(file = "baseline.svg", fix_text_size = FALSE)
baseline_gg
dev.off()


###
# get proportions of the topo types
tpi <- read_csv("src/data/baseline/stmNorthernForests tpi-fractions-table.csv", skip = 6) |>
    janitor::clean_names() |>
    filter(step == 0)

names(tpi)[3:5] <- c("prop_gully", "prop_slope", "prop_ridge")

tpi_lg <- tpi |>
    mutate(prop_gully = prop_gully / (256 ^ 2),
           prop_slope = prop_slope / (256 ^ 2),
           prop_ridge = prop_ridge / (256 ^ 2)) |>
    pivot_longer(3:5)

ggplot(data = tpi_lg) +
    geom_boxplot(aes(x = name, y = value))

