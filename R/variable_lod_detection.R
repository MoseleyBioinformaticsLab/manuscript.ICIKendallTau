create_variable_lod_samples = function(n_feature, n_sample,
                                       lod_values)
{
  # n_feature = 1000
  # n_sample = 100
  # lod_values = tibble::tibble(lod = 2.1, level = "basic")
  base_sample = rlnorm(n_feature, meanlog = 1, sdlog = 0.5)
  rep_data = add_uniform_noise(n_sample, base_sample, 0.2)
  colnames(rep_data) = paste0("S", seq_len(ncol(rep_data)))
  
  n_each_lod = ceiling(n_sample / nrow(lod_values))
  
  lod_vector = rep(lod_values$lod, each = n_each_lod)
  lod_vector = lod_vector[seq_len(n_sample)]
  
  lod_sample_tbl = tibble::tibble(lod = lod_vector, sample = colnames(rep_data))
  lod_sample_tbl = dplyr::left_join(lod_sample_tbl, lod_values, by = "lod")
  
  lod_matrix = matrix(lod_vector, nrow = nrow(rep_data), ncol = ncol(rep_data), byrow = TRUE)
  
  rep_lod = rep_data
  rep_lod[rep_lod < lod_matrix] = NA
  
  return(list(sample_data = rep_data,
              sample_lod = rep_lod,
              lod = lod_sample_tbl))
}

variable_lod_cor_everyway = function(sample_counts){
  # tmp_values = tar_read(var_lod_samples)
  # sample_counts = tmp_values$variable_cutoff
  
  min_value = min(sample_counts, na.rm = TRUE)
  impute_value = min_value / 2
  
  sample_counts_min = sample_counts
  sample_counts_min[is.na(sample_counts)] = impute_value
  
  ici_cor = ici_kendalltau(t(sample_counts), global_na = NA, scale_max = FALSE)$cor
  kt_na = kt_fast(sample_counts, use = "pairwise.complete.obs", return_matrix = TRUE)$tau
  kt_min = kt_fast(sample_counts_min, use = "pairwise.complete.obs", return_matrix = TRUE)$tau
  pearson_na = cor(sample_counts, method = "pearson", use = "pairwise.complete.obs")
  pearson_min = cor(sample_counts_min, method = "pearson", use = "pairwise.complete.obs")
  
  cor_vals = list(icikt = ici_cor,
                  kt_na = kt_na,
                  kt_min = kt_min,
                  pearson_na = pearson_na,
                  pearson_min = pearson_min
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

calculate_var_lod_correlation_diffs = function(var_lod_correlations)
{
  # tar_load(var_lod_correlations)
  reference_cor = var_lod_correlations$no_cutoff
  
  do_diffs = c("single_cutoff", "variable_cutoff")
  diff_cor = purrr::map(do_diffs, \(diff_id){
    # diff_id = do_diffs[1]
    purrr::map(names(reference_cor), \(cor_id){
      # cor_id = "icikt_na"
      diff_df = var_lod_each_cor(reference_cor[[cor_id]], var_lod_correlations[[diff_id]][[cor_id]])
      
      diff_df = diff_df |>
        dplyr::mutate(cor_method = cor_id,
                      lod_method = diff_id)
      diff_df
    }) |>
      purrr::list_rbind()
  }) |>
    purrr::list_rbind()
  
  diff_cor_lod = add_lod_info(diff_cor, var_lod_correlations$samples$cutoffs)
  
  diff_cor_lod
}

add_lod_info = function(diff_cor, cutoffs)
{
  # cutoffs = var_lod_correlations$samples$cutoffs
  
  cutoff_df = tibble::tibble(s1 = names(cutoffs), s2 = names(cutoffs), cutoff = cutoffs)
  uniq_cutoffs = tibble::tibble(cutoff = unique(cutoffs), level = c("low", "med", "high", "vhigh")) 
  cutoff_df = dplyr::left_join(cutoff_df, uniq_cutoffs, by = "cutoff")
  
  diff_cor_lod = dplyr::left_join(diff_cor, cutoff_df |> dplyr::transmute(s1 = s1, s1_cutoff = cutoff, s1_level = level), by = "s1")
  diff_cor_lod = dplyr::left_join(diff_cor_lod, cutoff_df |> dplyr::transmute(s2 = s2, s2_cutoff = cutoff, s2_level = level), by = "s2")
  
  diff_cor_lod = diff_cor_lod |>
    dplyr::mutate(compare_levels = glue::glue("{s1_level}-{s2_level}"))
  diff_cor_lod
}

var_lod_each_cor = function(ref_cor, in_cor)
{
  # ref_cor = reference_cor[["icikt_na"]]
  # in_cor = var_lod_correlations[["single_cutoff"]][["icikt_na"]]
  ref_long = cor_matrix_2_long_df(ref_cor) |>
    dplyr::rowwise() |>
    dplyr::mutate(comparison = paste0(sort(c(s1, s2)), collapse = ".")) |>
    dplyr::ungroup() |>
    dplyr::filter(!(s1 == s2), !duplicated(comparison))
  in_long = cor_matrix_2_long_df(in_cor) |>
    dplyr::rowwise() |>
    dplyr::mutate(comparison = paste0(sort(c(s1, s2)), collapse = ".")) |>
    dplyr::ungroup() |>
    dplyr::filter(!(s1 == s2), !duplicated(comparison))
  
  compare_cor = dplyr::left_join(ref_long[, c("cor", "s1", "s2", "comparison")], in_long[, c("cor", "comparison")], suffix = c("_ref", "_lod"), by = "comparison")
  compare_cor = compare_cor |>
    dplyr::mutate(ref_v_lod = cor_ref - cor_lod)
  compare_cor
  
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
