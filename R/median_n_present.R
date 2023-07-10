group_study = function(counts_info, id)
{
  # counts_info = tar_read(yeast_counts_info)
  # counts_info = tar_read(nsclc_counts_info)
  counts_df = tibble::as_tibble(counts_info$counts)
  if (is.null(rownames(counts_info$counts))) {
    counts_df$feature = seq_len(nrow(counts_df))
  } else {
    counts_df$feature = rownames(counts_info$counts)
  }
  
  
  counts_long = counts_df |>
    tidyr::pivot_longer(cols = -feature, names_to = "sample", values_to = "count")
  counts_long = dplyr::left_join(counts_long, counts_info$info[, c("sample", "treatment")], by = "sample")
  counts_long = counts_long |>
    dplyr::filter(!is.na(treatment))
  
  group_medians = counts_long |>
    dplyr::group_by(treatment, feature) |>
    dplyr::summarise(median = median(count),
                     median_present = median(count[count > 0]),
                     n_present = sum(count > 0)) |>
    dplyr::ungroup()
  
  group_median_min = group_medians |>
    dplyr::group_by(n_present, treatment) |>
    dplyr::summarise(min_median = min(median_present)) |>
    dplyr::ungroup() |>
    dplyr::mutate(dataset = id)
  
  group_medians$dataset = id
  
  return(list(medians = group_medians,
              median_min = group_median_min))
  
}


correlate_medians_n_present = function(grouped_medians,
                                       max_quantile = 0.5)
{
  # grouped_medians = tar_read(group_medians_yeast)
  median_npresent = grouped_medians$medians |>
    dplyr::mutate(log_median = log10(median_present)) |>
    dplyr::filter(!is.na(median_present))
  quantile_cut = 10^quantile(median_npresent$log_median, probs = max_quantile, na.rm = TRUE)
  median_under_q = median_npresent |>
    dplyr::filter(median <= quantile_cut)
  median_treatment = split(median_under_q, median_under_q$treatment)
  cor_median = purrr::imap(median_treatment, function(x, y){
    out_cor = suppressWarnings(ICIKendallTau::ici_kt(x$n_present, x$median_present, perspective = "global"))
    tibble::tibble(cor = out_cor["tau"], pvalue = out_cor["pvalue"], treatment = y)
  }) |>
    dplyr::bind_rows()
  
  median_min = grouped_medians$median_min |>
    dplyr::filter(!is.na(min_median))
  min_treatment = split(median_min, median_min$treatment)
  cor_min = purrr::imap(min_treatment, function(x, y){
    out_cor = suppressWarnings(ICIKendallTau::ici_kt(x$n_present, x$min_median, perspective = "global"))
    tibble::tibble(cor = out_cor["tau"], pvalue = out_cor["pvalue"], treatment = y)
  }) |>
    dplyr::bind_rows()
  return(list(medians = list(medians = grouped_medians$medians,
                      correlation = cor_median),
       median_min = list(medians = grouped_medians$median_min,
                         correlation = cor_min),
       quantile = quantile_cut)

  )
}
