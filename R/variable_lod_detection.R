create_variable_lod_samples = function(n_feature, n_sample,
                                       base_value = 2.1,
                                       variables = c(0.25, 0.5, 1, 1.5))
{
  # n_feature = 1000
  # n_sample = 100
  # variables = c(0.25, 0.5, 1, 1.5)
  base_sample = rlnorm(n_feature, meanlog = 1, sdlog = 0.5)
  rep_data = add_uniform_noise(n_sample, base_sample, 0.2)
  colnames(rep_data) = paste0("S", seq_len(ncol(rep_data)))
  
  base_value_matrix = matrix(base_value, nrow = nrow(rep_data), ncol = ncol(rep_data), byrow = TRUE)
  
  rep_base = rep_data
  rep_base[rep_data < base_value_matrix] = NA
  
  variable_cutoffs = rep(variables * base_value, each = n_sample / length(variables))
  
  variable_value_matrix = matrix(variable_cutoffs, nrow = nrow(rep_data), ncol = ncol(rep_data), byrow = TRUE)
  
  rep_variable = rep_data
  rep_variable[rep_data < variable_value_matrix] = NA
  
  return(list(base_cutoff = exp(base_value),
              variables = variables,
              cutoffs = exp(variable_cutoffs),
              no_cutoff = exp(rep_data),
              single_cutoff = exp(rep_base),
              variable_cutoff = exp(rep_variable)))
}

variable_lod_cor_everyway = function(sample_counts){
  # tmp_values = tar_read(var_lod_samples)
  # sample_counts = tmp_values$variable_cutoff
  
  min_value = min(sample_counts, na.rm = TRUE)
  impute_value = min_value / 2
  
  sample_counts_min = sample_counts
  sample_counts_min[is.na(sample_counts)] = impute_value
  
  ici_cor = ici_kendalltau(t(sample_counts), global_na = NA, scale_max = FALSE)$cor
  ici_min = ici_kendalltau(t(sample_counts_min), global_na = NA, scale_max = FALSE)$cor
  kt = kt_fast(sample_counts, use = "pairwise.complete.obs", return_matrix = TRUE)$tau
  kt_min = kt_fast(sample_counts_min, use = "pairwise.complete.obs", return_matrix = TRUE)$tau
  # this one should match the Gierlinski paper values for median correlations
  pearson_base_na = cor(sample_counts, method = "pearson", use = "pairwise.complete.obs")
  pearson_base_min = cor(sample_counts_min, method = "pearson", use = "pairwise.complete.obs")
  pearson_log1p = cor(log1p(sample_counts), method = "pearson", use = "pairwise.complete.obs")
  pearson_log1p_min = cor(log1p(sample_counts_min), method = "pearson", use = "pairwise.complete.obs")
  
  cor_vals = list(icikt_na = ici_cor,
                  icikt_min = ici_min,
                  kt_na = kt,
                  kt_min = kt_min,
                  pearson_base_na = pearson_base_na,
                  pearson_base_min = pearson_base_min,
                  pearson_log1p = pearson_log1p,
                  pearson_log1p_min = pearson_log1p_min
  )
  cor_vals
}



calculate_variable_correlations = function(var_lod_samples)
{
  # tar_load(var_lod_samples)
  
  nocutoff_correlations = variable_lod_cor_everyway(var_lod_samples$no_cutoff)
  singlecutoff_correlations = variable_lod_cor_everyway(var_lod_samples$single_cutoff)
  variablecutoff_correlations = variable_lod_cor_everyway(var_lod_samples$variable_cutoff)
  
  list(samples = var_lod_samples,
       no_cutoff = nocutoff_correlations,
       single_cutoff = singlecutoff_correlations,
       variable_cutoff = variablecutoff_correlations)
    
}
triple_check_scalemax_works = function(var_lod_samples)
{
  use_sample = var_lod_samples$single_cutoff
  
  ici_noscale = ici_kendalltau(t(use_sample), global_na = NA, scale_max = FALSE, diag_good = FALSE)
  
  testthat::expect_equal(ici_noscale$cor, ici_noscale$raw)
  
  ici_scale = ici_kendalltau(t(use_sample), global_na = NA, scale_max = TRUE, diag_good = FALSE)
  sum_lt = sum(ici_scale$cor <= ici_scale$raw)
  testthat::expect_equal(sum_lt, nrow(ici_scale$cor) * ncol(ici_scale$cor))
  NULL
}
