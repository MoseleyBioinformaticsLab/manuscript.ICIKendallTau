run_yeast_everyway = function(yeast_counts, yeast_completeness){
  
  ici_cor = ici_kendalltau(t(yeast_counts), global_na = c(NA, 0))$cor
  yeast_counts_na = yeast_counts
  yeast_counts_na[yeast_counts_na == 0] = NA
  pearson_base_nozero = cor(yeast_counts_na, method = "pearson", use = "pairwise.complete")
  pearson_base = cor(yeast_counts, method = "pearson", use = "pairwise.complete")
  pearson_log1p = cor(log1p(yeast_counts), method = "pearson", use = "pairwise.complete")
  log_counts = log(yeast_counts)
  log_counts[is.infinite(log_counts)] = NA
  pearson_log = cor(log_counts, method = "pearson", use = "pairwise.complete")
  
  cor_vals = list(icikt = ici_cor,
                  icikt_complete = ici_cor * yeast_completeness,
                  pearson_base = pearson_base,
                  pearson_base_nozero = pearson_base_nozero,
                  pearson_log1p = pearson_log1p,
                  pearson_log = pearson_log
                  )
  cor_vals
}

calculate_yeast_medians = function(yeast_cor, yeast_info){
  out_med = purrr::imap(yeast_cor, function(in_cor, cor_id){
    in_cor = in_cor[yeast_info$sample_rep, yeast_info$sample_rep]
    in_med = visualizationQualityControl::median_correlations(in_cor, yeast_info$sample)
    in_med$which = cor_id
    in_med
  })
}
