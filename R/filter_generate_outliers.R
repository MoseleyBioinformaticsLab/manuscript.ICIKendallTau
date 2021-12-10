filter_generate_outliers = function(counts, info, keep_num, sample_col, class_col){
  # counts = readRDS(here::here("data", "brainson_rnaseq201901_counts.rds"))
  # info = readRDS(here::here("data", "brainson_rnaseq201901_info.rds"))
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "tumor"
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  counts_completeness = pairwise_completeness(t(counts_filter))
  counts_cor = run_cor_everyway(counts_filter, counts_completeness)
  counts_medians = calculate_cor_medians(counts_cor, info[[sample_col]], info[[median_col]])
  
  counts_outliers = purrr::map(counts_medians, determine_outliers)
  counts_outliers$features = rownames(counts_filter)
  counts_outliers
}
