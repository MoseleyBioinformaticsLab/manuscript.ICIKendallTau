graph_median_min = function(median_correlation)
{
  # median_correlation = tar_read(correlate_medians_yeast)
  median_plots = create_median_plot(median_correlation$medians, median_correlation$quantile)
  min_plots = create_min_median_plot(median_correlation$median_min)
  list(median = median_plots,
       min = min_plots)
}

find_good_location = function(all_values, perc_range = 0.1)
{
  # range_values = c(-0.8652091, 18.1403589)
  # perc_range = 0.1
  
  range_values = range(all_values, na.rm = TRUE)
  range_diff = range_values[2] - range_values[1]
  
  use_loc = range_values[1] + (range_diff * perc_range)
  use_loc
}

create_median_plot = function(median_stuff, xlimit)
{
  # tmp = tar_read(correlate_medians_nsclc)
  # median_stuff = tmp$medians
  # xlimit = tmp$quantile
  n_treatments = length(unique(median_stuff$medians$treatment))
  if (n_treatments <= 3) {
    x_loc_fraction = 0.7
  } else {
    x_loc_fraction = 0.6
  }
  log_medians = median_stuff$medians |>
    dplyr::mutate(log_median = log10(median_present)) |>
    dplyr::filter(!is.na(log_median)) |>
    dplyr::filter(log_median <= log10(xlimit))
  use_locations = log_medians |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(cor_x = find_good_location(log_median, x_loc_fraction),
                     cor_y = find_good_location(n_present, 0.05))
  use_locations = dplyr::left_join(use_locations, median_stuff$correlation, by = "treatment")
  use_locations = use_locations |>
    dplyr::mutate(cor_label = paste0("\u03C4: ", format(cor, digits = 2)))
  
  nrow = 1
  if (n_treatments > 6) {
    nrow = 2
  }
  out_plot = log_medians |>
    ggplot(aes(x = log_median, y = n_present)) +
    geom_point() +
    geom_label(data = use_locations, aes(x = cor_x, y = cor_y, label = cor_label), hjust = 0) +
    facet_wrap(~ treatment, nrow = nrow, scales = "free") +
    theme(strip.background = NULL,
          strip.text.x = element_text(hjust = 0)) +
    labs(x = "Log10(Median-Present)", y = "N-Present")
  
  out_plot
}

create_min_median_plot = function(median_stuff)
{
  # tmp = tar_read(correlate_medians_yeast)
  # median_stuff = tmp$median_min
  # xlimit = tmp$quantile
  # 
  # tmp = tar_read(correlate_medians_typeandtumorculture)
  # median_stuff = tmp$median_min
  # xlimit = tmp$quantile
  n_treatments = length(unique(median_stuff$medians$treatment))
  if (n_treatments <= 3) {
    x_loc_fraction = 0.7
  } else {
    x_loc_fraction = 0.6
  }
  log_medians = median_stuff$medians |>
    dplyr::mutate(log_median = log10(min_median)) |>
    dplyr::filter(!is.na(log_median))
  use_locations = log_medians |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(cor_x = find_good_location(log_median, x_loc_fraction),
                     cor_y = find_good_location(n_present, 0.05))
  use_locations = dplyr::left_join(use_locations, median_stuff$correlation, by = "treatment")
  use_locations = use_locations |>
    dplyr::mutate(cor_label = paste0("\u03C4: ", format(cor, digits = 2)))
  
  nrow = 1
  if (n_treatments > 6) {
    nrow = 2
  }
  
  out_plot = log_medians |>
    ggplot(aes(x = log_median, y = n_present)) +
    geom_point() +
    geom_label(data = use_locations, aes(x = cor_x, y = cor_y, label = cor_label), hjust = 0) +
    facet_wrap(~ treatment, nrow = nrow, scales = "free") +
    theme(strip.background = NULL,
          strip.text.x = element_text(hjust = 0)) +
    labs(x = "Log10(Median-Present-Min)", y = "N-Present")
  
  out_plot
}
