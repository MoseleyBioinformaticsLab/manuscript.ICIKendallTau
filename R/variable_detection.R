create_variable_lod_samples = function(n_feature, n_sample)
{
  # n_feature = 1000
  # n_sample = 100
  base_sample = rlnorm(n_feature, meanlog = 1, sdlog = 0.5)
  rep_data = add_uniform_noise(n_sample, base_sample, 0.2)
  colnames(rep_data) = paste0("S", seq_len(ncol(rep_data)))
  actual_ici = ici_kendalltau(t(rep_data), scale_max = FALSE)
  
  exp_data = exp(rep_data)
  exp_ici = ici_kendalltau(t(exp_data), scale_max = FALSE)
  
  actual_pearson = cor(rep_data)
  exp_pearson = cor(exp_data)
  
  
  cutoff_dist = runif(n = n_sample, min = min(rep_data), max = median(rep_data))
  cutoff_matrix = matrix(cutoff_dist, nrow = nrow(rep_data), ncol = ncol(rep_data), byrow = TRUE)
  
  cutoff_rep = rep_data
  cutoff_rep[rep_data < cutoff_matrix] = NA
  
  cutoff_exp = exp(cutoff_rep)
  rep_exp = exp(rep_data)
  
  return(list(cutoff = cutoff_exp,
              full = rep_exp))
}
