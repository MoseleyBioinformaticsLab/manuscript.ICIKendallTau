create_large_replicate_samples = function(n_feature, n_sample)
{
  base_sample = rlnorm(n_feature, meanlog = 1, sdlog = 0.5)
  rep_data = add_uniform_noise(n_sample, base_sample, 0.2)
  colnames(rep_data) = paste0("S", seq_len(ncol(rep_data)))
  rep_data
}

create_variable_lod_samples = function(rep_data,
                                       lod_values, id)
{
  if (inherits(lod_values, "numeric")) {
    lod_values = tibble::tibble(lod = lod_values, level = id)
  } else if (!inherits(lod_values, "data.frame")) {
    stop("Please pass a vector or a data.frame")
  }
  
  n_lod = nrow(lod_values)
  
  n_sample = ncol(rep_data)
  n_each_lod = ceiling(n_sample / nrow(lod_values))
  
  lod_vector = rep(lod_values$lod, each = n_each_lod)
  lod_vector = lod_vector[seq_len(n_sample)]
  
  lod_sample_tbl = tibble::tibble(lod = lod_vector, sample = colnames(rep_data))
  lod_sample_tbl = dplyr::left_join(lod_sample_tbl, lod_values, by = "lod")
  
  lod_matrix = matrix(lod_vector, nrow = nrow(rep_data), ncol = ncol(rep_data), byrow = TRUE)
  
  rep_lod = rep_data
  rep_lod[rep_lod < lod_matrix] = NA
  
  rep_na = tibble::tibble(sample = colnames(rep_data), perc_na = apply(rep_lod, 2, \(x){sum(is.na(x)) / length(x)}))
  lod_sample_tbl = dplyr::left_join(lod_sample_tbl, rep_na, by = "sample")
  
  
  return(list(sample_data = rep_data,
              sample_lod = rep_lod,
              lod = lod_sample_tbl,
              n_lod = n_lod))
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
  # var_lod_samples = tar_read(vl_samples_med)
  
  ref_correlations = variable_lod_cor_everyway(var_lod_samples$sample_data)
  lod_correlations = variable_lod_cor_everyway(var_lod_samples$sample_lod)
  
  list(samples = var_lod_samples,
       lod = var_lod_samples$lod,
       n_lod = var_lod_samples$n_lod,
       reference_cor = ref_correlations,
       lod_cor = lod_correlations)
    
}

calculate_var_lod_correlation_diffs = function(var_lod_correlations)
{
  # var_lod_correlations = tar_read(vl_cor_med)
  reference_cor = var_lod_correlations$reference
  lod_cor = var_lod_correlations$lod_cor
  
  diff_cor = purrr::map(names(reference_cor), \(cor_id){
    diff_df = var_lod_each_cor(reference_cor[[cor_id]], lod_cor[[cor_id]])
      
      diff_df = diff_df |>
        dplyr::mutate(cor_method = cor_id)
      diff_df
  }) |>
    purrr::list_rbind()
  
  diff_cor_lod = add_lod_info(diff_cor, var_lod_correlations$lod)
  diff_cor_lod$n_lod = var_lod_correlations$n_lod
  
  diff_cor_lod
}

add_lod_info = function(diff_cor, lod_df)
{
  # lod_df = var_lod_correlations$lod
  
  
  diff_cor_lod = dplyr::left_join(diff_cor, lod_df |> dplyr::transmute(s1 = sample, s1_lod = lod, s1_level = level, s1_perc_na = perc_na), by = "s1")
  diff_cor_lod = dplyr::left_join(diff_cor_lod, lod_df |> dplyr::transmute(s2 = sample, s2_lod = lod, s2_level = level, s2_perc_na = perc_na), by = "s2")
  
  diff_cor_lod = diff_cor_lod |>
    dplyr::mutate(compare_levels = glue::glue("{s1_level}-{s2_level}"))
  diff_cor_lod
}

var_lod_each_cor = function(ref_cor, in_cor)
{
  # ref_cor = reference_cor[["icikt_na"]]
  # in_cor = var_lod_correlations[["single_cutoff"]][["icikt_na"]]
  ref_long = lod_cor_matrix_2_df(ref_cor) |>
    dplyr::mutate(comparison = glue::glue("{s1}.{s2}")) |>
    dplyr::filter(!(s1 == s2), !duplicated(comparison))
  in_long = lod_cor_matrix_2_df(in_cor) |>
    dplyr::mutate(comparison = glue::glue("{s1}.{s2}")) |>
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

lod_cor_matrix_2_df = function(in_matrix)
{
  comparisons = combn(colnames(in_matrix), 2)
  
  cor_vals = vector(mode = "double", length = ncol(comparisons))
  
  for (icol in seq_len(ncol(comparisons))) {
    cor_vals[icol] = in_matrix[comparisons[1, icol], comparisons[2, icol]]
  }
  
  cor_df = data.frame(s1 = comparisons[1, ], s2 = comparisons[2, ], cor = cor_vals)
  cor_df
}
