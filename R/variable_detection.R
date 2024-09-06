create_variable_lod_samples = function(n_feature, n_sample,
                                       variables = c(0.25, 0.5, 1, 1.5))
{
  # n_feature = 1000
  # n_sample = 100
  # variables = c(0.25, 0.5, 1, 1.5)
  base_sample = rlnorm(n_feature, meanlog = 1, sdlog = 0.5)
  rep_data = add_uniform_noise(n_sample, base_sample, 0.2)
  colnames(rep_data) = paste0("S", seq_len(ncol(rep_data)))
  
  base_value = 2.1
  
  base_value_matrix = matrix(base_value, nrow = nrow(rep_data), ncol = ncol(rep_data), byrow = TRUE)
  
  rep_base = rep_data
  rep_base[rep_data < base_value_matrix] = NA
  
  variable_cutoffs = rep(variables * base_value, each = n_sample / length(variables))
  
  variable_value_matrix = matrix(variable_cutoffs, nrow = nrow(rep_data), ncol = ncol(rep_data), byrow = TRUE)
  
  rep_variable = rep_data
  rep_variable[rep_data < variable_value_matrix] = NA
  
  return(list(cutoffs = exp(variable_cutoffs),
              no_cutoff = exp(rep_data),
              single_cutoff = exp(rep_base),
              variable_cutoff = exp(rep_variable)))
}
