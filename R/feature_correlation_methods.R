ici = function(counts_info, id, keep_num, sample_col, class_col)
{
  # counts_info = tar_read(yeast_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  tmp_out = ici_kendalltau(counts_filter, global_na = c(NA, 0))$cor
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "ici",
       full_id = id)
}

ici_completeness = function(counts_info, id, keep_num, sample_col, class_col)
{
  # counts_info = tar_read(yeast_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  counts_completeness = pairwise_completeness(counts_filter)
  tmp_out = ici_kendalltau(counts_filter, global_na = c(NA, 0))$cor * counts_completeness
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "ici_completeness",
       full_id = id)
}


kt = function(counts_info, id, keep_num, sample_col, class_col)
{
  # counts_info = tar_read(yeast_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  tmp_out = ici_kendalltau(counts_filter, global_na = c(NA), scale_max = FALSE, diag_good = FALSE)$cor
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "kt",
       full_id = id)
}


pearson_base_nozero = function(counts_info, id, keep_num, sample_col, class_col)
{
  
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  counts_filter_na = counts_filter
  counts_filter_na[counts_filter == 0] = NA
  # this one should match the Gierlinski paper values for median correlations
  tmp_out = cor(t(counts_filter_na), method = "pearson", use = "pairwise.complete")
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "pearson_base_nozero",
       full_id = id)
}

pearson_base = function(counts_info, id, keep_num, sample_col, class_col)
{
  
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  tmp_out = cor(t(counts_filter), method = "pearson", use = "pairwise.complete")
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "pearson_base",
       full_id = id)
}

pearson_log1p = function(counts_info, id, keep_num, sample_col, class_col)
{
  # counts_info = tar_read(yeast_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  
  
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  tmp_out = cor(log1p(t(counts_filter)), method = "pearson", use = "pairwise.complete")
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "pearson_log1p",
       full_id = id)
}

pearson_log = function(counts_info, id, keep_num, sample_col, class_col)
{
  
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  log_counts = log(counts_filter)
  log_counts[is.infinite(log_counts)] = NA
  tmp_out = cor(t(log_counts), method = "pearson", use = "pairwise.complete")
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "pearson_log",
       full_id = id)
}

