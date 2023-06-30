group_yeast_study = function(yeast_counts_info)
{
  #tar_load(yeast_counts_info)
  
  suppressMessages(library(DESeq2))
  yeast_counts = yeast_counts_info$counts
  
  yeast_des = DESeqDataSetFromMatrix(yeast_counts_info$counts, yeast_counts_info$info, design = ~ sample)
  yeast_des = estimateSizeFactors(yeast_des)
  yeast_norm = counts(yeast_des, normalized = TRUE)
  
  yeast_counts_df = tibble::as_tibble(yeast_norm)
  yeast_counts_df$feature = rownames(yeast_norm)
  
  yeast_counts_long = yeast_counts_df |>
    tidyr::pivot_longer(cols = -feature, names_to = "sample", values_to = "count")
  yeast_counts_long = dplyr::left_join(yeast_counts_long, yeast_counts_info$info, by = c("sample" = "sample_rep"))
  
  group_medians = yeast_counts_long |>
    dplyr::group_by(sample.y, feature) |>
    dplyr::summarise(median = median(count),
                     median_present = median(count[count > 0]),
                     n_present = sum(count > 0)) |>
    dplyr::ungroup() |>
    dplyr::mutate(group = sample.y) |>
    dplyr::select(-sample.y)
  
  group_median_min = group_medians |>
    dplyr::group_by(n_present, group) |>
    dplyr::summarise(min_median = min(median_present)) |>
    dplyr::ungroup()
  
  return(list(medians = group_medians,
              median_min = group_median_min))
  
}

group_brainsonrnaseq_study = function(brainsonrnaseq_counts,
                                      brainsonrnaseq_info)
{
  tar_load(brainsonrnaseq_counts)
  tar_load(brainsonrnaseq_info)
  suppressMessages(library(DESeq2))
  
  brainsonrnaseq_des = DESeqDataSetFromMatrix(round(brainsonrnaseq_counts), brainsonrnaseq_info, design = ~ tumor + type) 
  brainsonrnaseq_des = estimateSizeFactors(brainsonrnaseq_des)
  brainsonrnaseq_norm = counts(brainsonrnaseq_des, normalized = TRUE)
  
  brainsonrnaseq_norm_df = tibble::as_tibble(brainsonrnaseq_norm) |>
    dplyr::mutate(feature = rownames(brainsonrnaseq_norm))
  
  brainsonrnaseq_norm_long = brainsonrnaseq_norm_df |>
    tidyr::pivot_longer(-feature,
                        values_to = "count",
                        names_to = "sample")
  
  brainsonrnaseq_norm_long = dplyr::left_join(brainsonrnaseq_norm_long,
                                              brainsonrnaseq_info, by = "sample")
  
  group_medians = brainsonrnaseq_norm_long |>
    dplyr::mutate(group = paste0(type, ":", tumor)) |>
    dplyr::group_by(group, feature) |>
    dplyr::summarise(median = median(count),
                     median_present = median(count[count > 0]),
                     n_present = sum(count > 0)) |>
    dplyr::ungroup()
  
  group_median_min = group_medians |>
    dplyr::group_by(n_present, group) |>
    dplyr::summarise(min_median = min(median_present)) |>
    dplyr::ungroup()
  
  return(list(medians = group_medians,
              median_min = group_median_min))
  
}

group_nsclc_study = function(nsclc_peaks,
                             nsclc_medians,
                             nsclc_info)
{
  # tar_load(nsclc_peaks)
  # tar_load(nsclc_medians)
  # tar_load(nsclc_info)
  
  keep_samples = base::intersect(colnames(nsclc_peaks), nsclc_medians$sample)
  keep_samples = base::intersect(keep_samples, nsclc_info$sample)
  
  nsclc_info = nsclc_info |>
    dplyr::filter(sample %in% keep_samples) |>
    dplyr::mutate(disease_instrument = paste0(disease, ":", instrument))
  split_type_instrument = split(nsclc_info, nsclc_info$disease_instrument)
  
  nsclc_medians = nsclc_medians |>
    dplyr::filter(sample %in% keep_samples)
  
  nsclc_peaks = nsclc_peaks[, keep_samples]
  
  nsclc_norm_list = purrr::map(colnames(nsclc_peaks), function(.x){
    nsclc_peaks[, .x] / (dplyr::filter(nsclc_medians, sample %in% .x) |> dplyr::pull(median))
  })
  
  nsclc_norm = do.call(cbind, nsclc_norm_list)
  colnames(nsclc_norm) = colnames(nsclc_peaks)
  nsclc_norm[is.na(nsclc_norm)] = 0
  nsclc_norm_df = as_tibble(nsclc_norm)
  nsclc_norm_df$feature = seq_len(nrow(nsclc_norm_df))
  
  nsclc_norm_long = nsclc_norm_df |>
    tidyr::pivot_longer(-feature,
                        names_to = "sample",
                        values_to = "intensity")
  nsclc_norm_long = dplyr::left_join(nsclc_norm_long, nsclc_info[, c("sample", "instrument", "disease")])
  
  group_medians = nsclc_norm_long |>
    dplyr::mutate(group = paste(instrument, ":", disease)) |>
    dplyr::group_by(group, feature) |>
    dplyr::summarise(median_present = median(intensity[intensity > 0]),
                     n_present = sum(intensity > 0)) |>
    dplyr::ungroup()
  
  group_median_min = group_medians |>
    dplyr::group_by(n_present, group) |>
    dplyr::summarise(min_median = min(median_present)) |>
    dplyr::ungroup()
  
  return(list(medians = group_medians,
              median_min = group_median_min))
  
}

group_adenocarcinoma_study = function(adeno_data,
                                      adeno_info)
{
  # tar_load(adeno_data)
  # tar_load(adeno_info)
  suppressMessages(library(DESeq2))
  
  adeno_des = DESeqDataSetFromMatrix(adeno_data, adeno_info, design = ~ tissue_type) 
  adeno_des = estimateSizeFactors(adeno_des)
  adeno_norm = counts(adeno_des, normalized = TRUE)
  
  adeno_counts_df = tibble::as_tibble(adeno_norm)
  adeno_counts_df$feature = rownames(adeno_norm)
  
  adeno_counts_long = adeno_counts_df |>
    tidyr::pivot_longer(cols = -feature, names_to = "sample_id2", values_to = "count")
  adeno_counts_long = dplyr::left_join(adeno_counts_long, adeno_info[, c("sample_id2", "tissue_type")], by = c("sample_id2" = "sample_id2"))
  
  group_medians = adeno_counts_long |>
    dplyr::group_by(tissue_type, feature) |>
    dplyr::summarise(median = median(count),
                     median_present = median(count[count > 0]),
                     n_present = sum(count > 0)) |>
    dplyr::ungroup() |>
    dplyr::mutate(group = tissue_type) |>
    dplyr::select(-tissue_type)
  
  group_median_min = group_medians |>
    dplyr::group_by(n_present, group) |>
    dplyr::summarise(min_median = min(median_present)) |>
    dplyr::ungroup()
  
  return(list(medians = group_medians,
              median_min = group_median_min))
  
  
}

group_mwtab_study = function(mwtab_normalized)
{
  # tar_load(mwtab_normalized)
  norm_data = mwtab_normalized$normalized_data
  sample_info = mwtab_normalized$sample_info
  
  use_factors = grep("^factor", names(sample_info), value = TRUE)
  group_factor = split(sample_info$sample_id, sample_info[, use_factors])
  group_factor = purrr::map(group_factor, function(.x){
    .x[.x %in% colnames(norm_data)]
  })
  
  group_data = purrr::imap(group_factor, function(in_group, group_id){
    use_samples = c(in_group, "metabolite")
    group_normalized = norm_data |>
      dplyr::select(tidyselect::all_of(use_samples)) |>
      dplyr::mutate(group = group_id)
    group_normalized
  }) |>
    dplyr::bind_rows()
  
  group_long = group_data |>
    tidyr::pivot_longer(c(-metabolite, -group),
                        names_to = "sample_id",
                        values_to = "intensity")
  
  group_long$intensity[is.na(group_long$intensity)] = 0
  group_medians = group_long |>
    dplyr::group_by(group, metabolite) |>
    dplyr::summarise(n_present = sum(intensity > 0),
                     median_present = median(intensity[intensity > 0]),
                     median = median(intensity),
                     mean = mean(intensity),
                     mean_present = mean(intensity[intensity > 0])) |>
    dplyr::ungroup()
  group_median_min = group_medians |>
    dplyr::group_by(group, n_present) |>
    dplyr::summarise(min_median = min(median_present)) |>
    dplyr::ungroup()
  
  return(list(medians = group_medians,
              median_min = group_median_min))
}

correlate_medians_n_present = function(grouped_medians,
                                       max_quantile = 0.5)
{
  # grouped_medians = tar_read(nsclc_grouped)
  median_npresent = grouped_medians$medians |>
    dplyr::mutate(log_median = log2(median_present)) |>
    dplyr::filter(!is.na(median_present))
  quantile_cut = quantile(median_npresent$log_median, probs = max_quantile, na.rm = TRUE)
  median_under_q = median_npresent |>
    dplyr::filter(log_median <= quantile_cut)
  median_group = split(median_under_q, median_under_q$group)
  cor_median = purrr::imap(median_group, function(x, y){
    out_cor = suppressWarnings(ICIKendallTau::ici_kt(x$n_present, x$median_present, perspective = "global"))
    tibble::tibble(cor = out_cor["tau"], pvalue = out_cor["pvalue"], group = y)
  }) |>
    dplyr::bind_rows()
  
  median_min = grouped_medians$median_min |>
    dplyr::filter(!is.na(min_median))
  min_group = split(median_min, median_min$group)
  cor_min = purrr::imap(min_group, function(x, y){
    out_cor = suppressWarnings(ICIKendallTau::ici_kt(x$n_present, x$min_median, perspective = "global"))
    tibble::tibble(cor = out_cor["tau"], pvalue = out_cor["pvalue"], group = y)
  }) |>
    dplyr::bind_rows()
  list(medians = cor_median,
       median_min = cor_min)
}
