## Load your packages, e.g. library(targets).
tar_source("./packages.R")

## Load your R files
tar_source(files = "R")

## tar_plan supports drake-style targets and also tar_target()
tar_plan(

  # small theoretical example --------
  x = seq(1, 10),
  y = seq(1, 10),
  y2 = seq(10, 1),
  where_na = create_na_indices(20),
  
  all_kt = all_kendalltau(x, y, y2, where_na),
  
  positive_kt = compare_positive_kt(x, y, where_na),
  negative_kt = compare_negative_kt(x, y2, where_na),
  positive_pearson = compare_positive_pearson(x, y, where_na),
  negative_pearson = compare_negative_pearson(x, y2, where_na),
  
  positive_kendall = compare_positive_pearson(x, y, where_na, method = "kendall"),
  negative_kendall = compare_negative_pearson(x, y2, where_na, method = "kendall"),
  
  # bigger, more realistic example ---------
  realistic_sample_1 = create_sample(n = 1000),
  realistic_sample_2 = create_sample(n = 1000),
  realistic_neg_sample = sort(realistic_sample_2, decreasing = TRUE),
  realistic_na = create_random_na(),
  realistic_positive_kt = compare_positive_kt(realistic_sample_1, realistic_sample_2, realistic_na),
  realistic_negative_kt = compare_negative_kt(realistic_sample_1, realistic_neg_sample, realistic_na),
  
  realistic_positive_pearson = compare_positive_pearson(realistic_sample_1, realistic_sample_2, realistic_na),
  realistic_negative_pearson = compare_negative_pearson(realistic_sample_1, realistic_neg_sample, realistic_na),
  
  realistic_positive_kendall = compare_positive_kendall(realistic_sample_1, realistic_sample_2, realistic_na),
  realistic_negative_kendall = compare_negative_kendall(realistic_sample_1, realistic_neg_sample, realistic_na),
  
  left_censored_samples = create_lc_samples(),
  left_censored_cor = left_censor_correlate(left_censored_samples),
  random_censored_cor = random_censor_correlate(left_censored_samples),
  logtransform_censored_cor = lt_left_censor_correlate(left_censored_samples),

  subsample = sample(1000, 50),
  left_sample_cor = left_censor_correlate(left_censored_samples[subsample, ]),
  random_sample_cor = random_censor_correlate(left_censored_samples[subsample, ], n_na = seq(0, 12, 2)),
  logtransform_sample_cor = lt_left_censor_correlate(left_censored_samples[subsample, ]),

  # loading real data ---------
  tar_target(adenocarcinoma_file,
             here::here("data/recount_adenocarcinoma_count_info.rds"),
             format = "file"),
  adenocarcinoma_data_info = readRDS(adenocarcinoma_file),
  
  tar_target(barton_yeast_file,
             here::here("data/barton_yeast_counts_info.rds"),
             format = "file"),
  yeast_data_info = readRDS(barton_yeast_file),
  
  tar_target(brainson_tumorcultures_file,
             here::here("data/brainsonrnaseq_type_counts_info.rds"),
             format = "file"),
  tumorcultures_data_info = readRDS(brainson_tumorcultures_file),
  
  tar_target(mwtab_micediabetics_file,
             here::here("data/mwtab_st000017_an000034_count_info.rds"),
             format = "file"),
  micediabetics_data_info = readRDS(mwtab_micediabetics_file),
  
  tar_target(rcsirm_nsclc_file,
             here::here("data/nsclc_count_info.rds"),
             format = "file"),
  nsclc_data_info = readRDS(rcsirm_nsclc_file),
  
  
  # running a single core to measure performance aspects ----------
  single_core_perf = run_single_cor(),
  complexity_figure = create_complexity_figure(single_core_perf),
  
  tumorcultures_n = c(1, 0.25, 0.5, 0.75, 0.99),
  tar_target(tumorcultures_outliers,
             filter_generate_outliers(tumorcultures_counts, tumorcultures_info, brainson_n, "sample", "tumor"),
             pattern = map(brainson_n),
             iteration = "list"),
  tumorcultures_single = get_single_outlier(tumorcultures_outliers),

  tar_target(tumorcultures_outliers_alt,
             filter_generate_outliers(tumorcultures_counts, tumorcultures_info, brainson_n, "sample", c("type", "tumor")),
             pattern = map(brainson_n),
             iteration = "list"),
  
  yeast_paper_outliers = c("WT.21", "WT.22", "WT.25", "WT.28", "WT.34", "WT.36",
                           "Snf2.06", "Snf2.13", "Snf2.25", "Snf2.35"),
  yeast_counts_info = readRDS(here::here("data", "yeast_counts_info.rds")),
  
  yeast_n = c(1, 0.25, 0.5),
  tar_target(yeast_outliers,
             filter_generate_outliers(yeast_counts_info$counts, yeast_counts_info$info, yeast_n, "sample_rep", "sample"),
             pattern = map(yeast_n),
             iteration = "list"),
  yeast_single = get_single_outlier(yeast_outliers),
  
  adeno_info = readRDS(here::here("data", "transcript_info.rds")),
  adeno_data = transcript_data,
  
  tar_target(adeno_outliers,
             filter_generate_outliers(adeno_data, adeno_info, yeast_n, "sample_id2", "tissue_type"),
             pattern = map(yeast_n),
             iteration = "list"),
  adeno_single = get_single_outlier(adeno_outliers),
  
  multi_samples = seq(5, 95, 5),
  tar_target(select_ss_multi,
             select_samples(yeast_counts_info$counts, multi_samples),
             pattern = map(multi_samples),
             iteration = "list"
  ),
  
  tar_target(run_multi,
             run_big_samples(select_ss_multi),
             pattern = map(select_ss_multi),
             iteration = "list"),
  
  tar_target(time_multi,
             get_run_time(run_multi),
             pattern = map(run_multi)),
  
  
  mwtab_grouped_correlation = correlate_medians_n_present(mwtab_grouped),
  
  nsclc_grouped = group_nsclc_study(nsclc_peaks,
                                    nsclc_medians,
                                    nsclc_info),
  nsclc_grouped_correlation = correlate_medians_n_present(nsclc_grouped),
  
  adeno_grouped = group_adenocarcinoma_study(adeno_data,
                                             adeno_info),
  adeno_grouped_correlation = correlate_medians_n_present(adeno_grouped),
  yeast_grouped = group_yeast_study(yeast_counts_info),
  yeast_grouped_correlation = correlate_medians_n_present(yeast_grouped),
  brainsonrnaseq_grouped = group_brainsonrnaseq_study(brainsonrnaseq_counts,
                                                      brainsonrnaseq_info),
  brainsonrnaseq_grouped_correlation = correlate_medians_n_present(brainsonrnaseq_grouped),
  
  
  tar_render(supp_materials,
              "doc/supplemental_materials.Rmd"),
  tar_render(supp_tables,
             "doc/supplemental_tables.Rmd"),
  tar_render(manuscript,
           "doc/ici_kt_manuscript.Rmd")

  
# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style

)
