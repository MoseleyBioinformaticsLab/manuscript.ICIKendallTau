filter_generate_outliers = function(counts, info, keep_num, sample_col, class_col){
  # counts_info = tar_read(yeast_counts_info)
  # counts = counts_info$counts
  # info = counts_info$info
  # 
  # counts_info = tar_read(nsclc_counts_info)
  # counts = counts_info$counts
  # info = counts_info$info
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
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
  
  counts_outliers = purrr::map_df(counts_medians, determine_outliers)
  counts_outliers$keep_num = keep_num
  list(outliers = counts_outliers,
       features = rownames(counts))
}


get_single_outlier = function(outlier_list)
{
  # outlier_list = tar_read(yeast_outliers)
  n_frac = purrr::map_dbl(outlier_list, \(x){x$outliers$keep_num[1]})
  outlier_list[[which(n_frac == 1)]]$outliers
}

rbind_outliers = function(outlier_list)
{
  all_outliers = purrr::map(outlier_list, \(x){x$outliers}) |>
    dplyr::bind_rows()
  all_outliers
}
