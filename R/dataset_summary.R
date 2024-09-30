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
  info_by_treatment = info |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(n_rep = dplyr::n())
  n_treatment = nrow(info_by_treatment)
  n_rep_treatment = paste0(info_by_treatment$n_rep, collapse = ", ")
  
  tibble::tibble(Features = n_feature,
                 Samples = n_sample,
                 Treatments = n_treatment,
                 Replicates = n_rep_treatment)
}

create_lc_summary = function(lc_test_result)
{
  # lc_test_result = tar_read(lc_test_ratstamina)
  missingness = lc_test_result$missingness |>
    dplyr::filter((treatment %in% "All")) |>
    dplyr::mutate(perc_value = format(perc_na, digits = 2))
  n_perc_condition = paste0(missingness$perc_value, collapse = ", ")
  
  lc_res = lc_test_result$stats |>
    dplyr::transmute(Dataset = id,
                     `Total` = missingness$n,
                     `Missing` = missingness$n_na,
                     `% Missing` = missingness$perc_na,
                     `Trials` = parameter,
                     `Success` = statistic,
                     `Estimate` = estimate,
                     `CI` = paste0(format(conf.low, digits = 2, nsmall = 2), " - ", format(conf.high, digits = 2, nsmall = 2)),
                     `P-Value` = dplyr::case_when(
                       p.value < 2e-16 ~ 0,
                       TRUE ~ p.value
                     ))
    
  lc_res
}

format_exponent_flextable = function(in_table, i, j)
{
  actual_value = in_table$body$dataset[i, j]
  split_value = strsplit(actual_value, "e")[[1]]
  new_number = paste0(split_value[1], "×10")
  out_table = compose(in_table, i = i, j = j, part = "body",
                      value = as_paragraph(new_number, as_sup(split_value[2])))
  out_table
}

check_zeros = function(in_df, j)
{
  in_df[[j]] = gsub("^0.*", "0", in_df[[j]])
  in_df
}
