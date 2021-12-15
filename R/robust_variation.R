calculate_variation = function(median_cor_df){
  suppressMessages(median_cor_df %>%
    dplyr::group_by(method, sample_class, keep_num) %>%
    dplyr::summarise(median = median(med_cor), mad = mad(med_cor)))
}

calculate_differences = function(median_stats, groups = NULL){
  if (is.null(groups)) {
    groups = unique(median_stats$sample_class)
  }
  median_stats %>%
    dplyr::filter(sample_class %in% groups) %>%
    dplyr::group_by(method, keep_num) %>%
    dplyr::summarise(med_diff = median[1] - median[2], mad_diff = mad[1] - mad[2]) %>%
      tidyr::pivot_longer(cols = c(med_diff, mad_diff),
                          names_to = "which_diff",
                          values_to = "diff")
}
