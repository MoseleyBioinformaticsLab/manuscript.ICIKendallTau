run_cor_everyway = function(sample_counts, sample_completeness){
  
  ici_cor = ici_kendalltau(t(sample_counts), global_na = c(NA, 0))$cor
  ici_cor_kt = ici_kendalltau(t(sample_counts), global_na = c(NA), scale_max = FALSE, diag_good = FALSE)
  sample_counts_na = sample_counts
  sample_counts_na[sample_counts_na == 0] = NA
  pearson_base_nozero = cor(sample_counts_na, method = "pearson", use = "pairwise.complete")
  pearson_base = cor(sample_counts, method = "pearson", use = "pairwise.complete")
  pearson_log1p = cor(log1p(sample_counts), method = "pearson", use = "pairwise.complete")
  log_counts = log(sample_counts)
  log_counts[is.infinite(log_counts)] = NA
  pearson_log = cor(log_counts, method = "pearson", use = "pairwise.complete")
  
  cor_vals = list(icikt = ici_cor,
                  icikt_complete = ici_cor * sample_completeness,
                  pearson_base = pearson_base,
                  pearson_base_nozero = pearson_base_nozero,
                  pearson_log1p = pearson_log1p,
                  pearson_log = pearson_log,
                  kt_base = ici_cor_kt
                  )
  cor_vals
}

calculate_cor_medians = function(sample_cor, sample_ids, sample_classes){
  out_med = purrr::imap(sample_cor, function(in_cor, cor_id){
    in_cor = in_cor[sample_ids, sample_ids]
    in_med = visualizationQualityControl::median_correlations(in_cor, sample_classes)
    in_med$which = cor_id
    in_med
  })
}
