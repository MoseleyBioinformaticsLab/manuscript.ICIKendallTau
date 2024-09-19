create_lod_run_df = function(lod_combinations, lod_levels)
{
  run_df = purrr::map(lod_combinations, \(in_comb){
    # in_comb = lod_combinations[[3]]
    purrr::map(seq_len(ncol(in_comb)), \(in_col){
      use_comb = in_comb[, in_col]
      tibble::tibble(multipliers = list(lod_levels$lod[use_comb]),
                     id = paste0(lod_levels$level[use_comb], collapse = "_"))
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows()
  run_df
}

# an alternative method to limit of detection would be to use a varying dynamic
# range, where we vary from 3 - 2 orders of magnitude, and base it on the highest
# observed value. May be able to use a uniform random distribution.
# 
# Also consider the differences not just as an absolute difference, but as a relative
# difference, where we take the difference / true correlation.

create_large_replicate_samples = function(n_feature, n_sample)
{
  base_sample = rlnorm(n_feature, meanlog = 1, sdlog = 0.5)
  # should we increase meanlog = 3, sdlog = 1.5, then noise = 0.2
  rep_data = add_uniform_noise(n_sample, base_sample, 0.2)
  colnames(rep_data) = paste0("S", stringr::str_pad(seq_len(ncol(rep_data)), width = 3, pad = "0"))
  
  base_lod = quantile(rep_data, 0.3)
  names(base_lod) = NULL
  list(data = rep_data,
       lod = base_lod)
}

create_variable_lod_samples = function(rep_info,
                                       lod_values, id)
{
  #lod_values = lod_values[[1]]
  rep_data = rep_info$data
  base_lod = rep_info$lod
  n_lod = length(lod_values)
  n_each_lod = 100
  n_sample = n_lod * n_each_lod
  
  lod_vector = rep(lod_values * base_lod, each = n_each_lod)
  lod_vector = lod_vector[seq_len(n_sample)]
  lod_df = tibble::tibble(multipliers = lod_values, lod = lod_values * base_lod)
  
  sampleid_sample = sample(colnames(rep_data), n_sample, replace = FALSE)
  
  use_data = rep_data[, sampleid_sample]
  lod_sample_tbl = tibble::tibble(lod = lod_vector, sample = colnames(use_data))
  lod_sample_tbl = dplyr::left_join(lod_sample_tbl, lod_df, by = "lod")
  
  lod_matrix = matrix(lod_vector, nrow = nrow(use_data), ncol = ncol(use_data), byrow = TRUE)
  
  use_lod = use_data
  use_lod[use_data < lod_matrix] = NA
  
  use_na = tibble::tibble(sample = colnames(use_data), perc_na = apply(use_lod, 2, \(x){sum(is.na(x)) / length(x)}))
  lod_sample_tbl = dplyr::left_join(lod_sample_tbl, use_na, by = "sample")
  
  
  return(list(sample_data = use_data,
              sample_lod = use_lod,
              lod = lod_sample_tbl,
              n_lod = n_lod,
              lod_id = id))
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
  impute_value = min(var_lod_samples$sample_lod, na.rm = TRUE) / 2
  ref_correlations = variable_lod_cor_everyway(var_lod_samples$sample_data)
  lod_correlations = variable_lod_cor_everyway(var_lod_samples$sample_lod)
  
  list(samples = var_lod_samples,
       lod = var_lod_samples$lod,
       n_lod = var_lod_samples$n_lod,
       reference_cor = ref_correlations,
       lod_cor = lod_correlations,
       impute_value = impute_value)
    
}

calculate_var_lod_correlation_diffs = function(var_lod_correlations)
{
  # var_lod_correlations = tar_read(vl_cor_med)
  # var_lod_correlations = tar_read(vl_cor_0.5_1_1.25_1.5)
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
  diff_cor_lod$impute_value = var_lod_correlations$impute_value
  diff_cor_lod$lod_id = var_lod_correlations$samples$lod_id
  
  diff_cor_lod
}

add_lod_info = function(diff_cor, lod_df)
{
  # lod_df = var_lod_correlations$lod
  
  
  diff_cor_lod = dplyr::left_join(diff_cor, lod_df |> dplyr::transmute(s1 = sample, s1_lod = lod, s1_multiplier = multipliers, s1_perc_na = perc_na), by = "s1")
  diff_cor_lod = dplyr::left_join(diff_cor_lod, lod_df |> dplyr::transmute(s2 = sample, s2_lod = lod, s2_multiplier = multipliers, s2_perc_na = perc_na), by = "s2")
  
  diff_cor_lod = diff_cor_lod |>
    dplyr::mutate(compare_levels = glue::glue("{s1_multiplier}::{s2_multiplier}"))
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

lod_cor_sample_matrix_2_df = function(lod_cor)
{
  # lod_cor = tar_read(vl_cor_all)
  sample_matrix = lod_cor$samples$sample_data
  sample_df = sample_matrix_2_df(sample_matrix)
  sample_df$lod_id = lod_cor$samples$sample_id
  sample_df = dplyr::left_join(sample_df, lod_cor$lod, by = "sample")
  sample_df$imputed_value = lod_cor$impute_value
  sample_df$which = "full"
  
  lod_matrix = lod_cor$samples$sample_lod
  lod_df = sample_matrix_2_df(lod_matrix)
  lod_df$lod_id = lod_cor$samples$sample_id
  lod_df = dplyr::left_join(lod_df, lod_cor$lod, by = "sample")
  lod_df$imputed_value = lod_cor$impute_value
  lod_df$which = "lod"
  
  dplyr::bind_rows(sample_df, lod_df)
}

sample_matrix_2_df = function(sample_matrix)
{
  n_tot = nrow(sample_matrix) * ncol(sample_matrix)
  total_value = vector("numeric", n_tot)
  if (is.null(rownames(sample_matrix))) {
    row_names = paste0("f", seq_len(nrow(sample_matrix)))
  } else {
    row_names = rownames(sample_matrix)
  }
  
  row_id = vector("character", n_tot)
  col_id = vector("character", n_tot)
  i_entry = 1
  for (icol in colnames(sample_matrix)) {
    for (irow in seq_len(nrow(sample_matrix))) {
      row_id[i_entry] = row_names[irow]
      col_id[i_entry] = icol
      total_value[i_entry] = sample_matrix[irow, icol]
      i_entry = i_entry + 1
    }
  }
  tibble::tibble(feature = row_id, sample = col_id, value = total_value)
  
}

calculate_cor_diff_summaries = function(vl_cor_diff)
{
  # vl_cor_diff = tar_read(vl_cor_diff_0.5_1.5)
  vl_summary = vl_cor_diff |>
    dplyr::group_by(cor_method, compare_levels) |>
    dplyr::summarise(ref_v_lod_median = median(abs(ref_v_lod)),
                     s1_percna_median = median(s1_perc_na),
                     s2_percna_median = median(s2_perc_na),
                     n_lod = n_lod[1],
                     lod_id = lod_id[1])
  vl_summary
}
