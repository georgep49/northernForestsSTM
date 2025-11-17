# This takes the raw NL outpout and parses it to final state
# plus some summary stats
parse_to_state <- function(full_dat)
{
    full_dat |>
        filter(step == 0 | step == 300) |>
        mutate(abundances = str_remove_all(abundances, "\\[|\\]")) |>
        separate_wider_delim(abundances, delim = " ", names = class_names) |>
        mutate(text_data = str_remove_all(abundances_by_topo, "\\[|\\]")) |>  # remove brackets
        separate_wider_delim(cols = text_data, names = class_names_topo, delim = " ", too_few = "align_start") |>
        mutate(across(starts_with("prop_"), ~as.numeric(.x))) |>
    left_join(fireLHC_size, by = c("siminputrow" = "ensemble")) |>
    rowwise() |>
    mutate(
        prop_ksh = sum(prop_kshK, prop_kshNoK, prop_kshP),
        prop_yfor = sum(prop_yfK, prop_yfNok, prop_yfP),
        prop_ofor = sum(prop_old, prop_oldP),
        ) |>
    ungroup() |>
    group_by(siminputrow) |>   # change in forest proportion over run (log response ratio)
    mutate(delta_ofor = log((prop_ofor[step == 300] + 1) / prop_ofor[step == 0])) |>
    ungroup() |>
    filter(step == 300)

# +1 on the delta_ofor to deal with divde by zero errors
}