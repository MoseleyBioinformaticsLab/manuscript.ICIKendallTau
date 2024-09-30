calculate_max_nomax = function(counts_info)
{
  # counts_info = tar_read(yeast_counts_info)
  # counts = counts_info$counts
  # info = counts_info$info
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  # 
  counts = counts_info$counts
  info = counts_info$info
  keep_num = 1
  sample_col = "sample"
  class_col = "treatment"
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  
  ici_cor = ici_kendalltau(t(counts_filter), global_na = c(NA, 0))$cor
  ici_nomax = ici_kendalltau(t(counts_filter), global_na = c(NA, 0), scale_max = FALSE)$cor
  
  list(max = ici_cor,
       nomax = ici_nomax,
       info = info)
}


calculate_diffs = function(max_nomax)
{
  # tar_load(max_nomax)
  max_vals = max_nomax$max
  max_vals[upper.tri(max_vals, diag = TRUE)] = NA
  nomax_vals = max_nomax$nomax
  nomax_vals[upper.tri(nomax_vals, diag = TRUE)] = NA
  
  max_nomax_ratio = max_vals / nomax_vals
  uniq_vals = unique(as.vector(max_nomax_ratio))
  uniq_vals = uniq_vals[!is.na(uniq_vals)]
  uniq_str = format(uniq_vals, digits = 6)
  
  sample_ids = max_nomax$info$sample
  sample_classes = max_nomax$info$treatment
  names(sample_classes) = max_nomax$info$sample
  out_med = purrr::imap(max_nomax[c("max", "nomax")], function(in_cor, cor_id){
    intersect_ids = base::intersect(colnames(in_cor), sample_ids)
    keep_ids = sample_ids %in% intersect_ids
    use_ids = sample_ids[keep_ids]
    use_classes = sample_classes[keep_ids]
    
    in_cor = in_cor[use_ids, use_ids]
    in_med = visualizationQualityControl::median_correlations(in_cor, use_classes)
    in_med$which = cor_id
    in_med
  })
  
  all_med = out_med |> purrr::list_rbind()
  
  wide_med = all_med |>
    tidyr::pivot_wider(names_from = "which", values_from = "med_cor")
  
  max_nomax_plot = wide_med |>
    ggplot(aes(x = nomax, y = max)) +
    geom_abline(slope = uniq_vals, color = "red") +
    geom_point(size = 2) +
    coord_equal() +
    labs(x = "Un-Scaled Median Correlations", y = "Scaled Median Correlations",
         subtitle = glue::glue("Scaled / Un-Scaled Ratio: {uniq_str}"))
  max_nomax_plot
}
