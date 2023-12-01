create_dataset_summary = function(counts_info)
{
  # counts_info = tar_read(yeast_counts_info)
  counts = counts_info$counts
  info = counts_info$info
  
  keep_samples = base::intersect(colnames(counts), info$sample)
  counts = counts[, keep_samples]
  info = info |>
    dplyr::filter(sample %in% keep_samples)
  
  n_feature = nrow(counts)
  n_sample = ncol(counts)
  info_by_condition = info |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(n_rep = dplyr::n())
  n_condition = nrow(info_by_condition)
  n_rep_condition = paste0(info_by_condition$n_rep, collapse = ", ")
  
  tibble::tibble(Features = n_feature,
                 Samples = n_sample,
                 Conditions = n_condition,
                 Replicates = n_rep_condition)
}
