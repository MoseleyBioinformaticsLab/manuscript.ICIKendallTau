run_yeast_everyway = function(yeast_counts){
  ici_cor = ici_kendalltau(t(yeast_counts), global_na = c(NA, 0))$cor
  pearson_base = cor(yeast_counts, method = "pearson", use = "pairwise.complete")
  pearson_log1p = cor(log1p(yeast_counts), method = "pearson", use = "pairwise.complete")
  log_counts = log(yeast_counts)
  log_counts[is.infinite(log_counts)] = NA
  pearson_log = cor(log_counts, method = "pearson", use = "pairwise.complete")
  
  cor_vals = list(icikt = ici_cor,
                  pearson_base = pearson_base,
                  pearson_log1p = pearson_log1p,
                  pearson_log = pearson_log
                  )
  cor_vals
}
