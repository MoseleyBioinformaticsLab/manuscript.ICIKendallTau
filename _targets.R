## Load your packages, e.g. library(targets). -----
tar_source("./packages.R")

## Load your R files -------
tar_source(files = "R")

# creating mapping tibble --------
dataset_variables = tibble::tibble(character = c("yeast_counts_info",
                                                 "adenocarcinoma_counts_info",
                                  "egfrgenotype_counts_info",
                                  "egfrgenotypetumorculture_counts_info",
                                  "ratstamina_counts_info",
                                  "nsclc_counts_info")) |>
  dplyr::mutate(id = gsub("_.*", "", character),
                sym = rlang::syms(character))


small_realistic_examples = tar_plan(

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
  logtransform_sample_cor = lt_left_censor_correlate(left_censored_samples[subsample, ])
)

loading_real_data = tar_plan(
  # loading real data ---------
  tar_target(adenocarcinoma_file,
             here::here("data/recount_adenocarcinoma_count_info.rds"),
             format = "file"),
  adenocarcinoma_counts_info = readRDS(adenocarcinoma_file),
  
  tar_target(barton_yeast_file,
             here::here("data/barton_yeast_counts_info.rds"),
             format = "file"),
  yeast_counts_info = readRDS(barton_yeast_file),
  
  tar_target(brainson_egfrgenotype_file,
             here::here("data/brainsonrnaseq_type_counts_info.rds"),
             format = "file"),
  egfrgenotype_counts_info = readRDS(brainson_egfrgenotype_file),
  
  tar_target(brainson_tumorcultures_file,
             here::here("data/brainsonrnaseq_type_counts_info_tumor.rds"),
             format = "file"),
  egfrgenotypetumorculture_counts_info = readRDS(brainson_tumorcultures_file),
  
  tar_target(mwtab_ratstamina_file,
             here::here("data/mwtab_st000017_an000034_count_info.rds"),
             format = "file"),
  ratstamina_counts_info = readRDS(mwtab_ratstamina_file),
  
  tar_target(rcsirm_nsclc_file,
             here::here("data/nsclc_count_info.rds"),
             format = "file"),
  nsclc_counts_info = readRDS(rcsirm_nsclc_file)
)

# limit of detection examples --------
limit_of_detection_map = tar_map(dataset_variables,
                                 names = id,
                                 tar_target(group_medians,
                                            group_study(sym, id)),
                                 tar_target(correlate_medians,
                                            correlate_medians_n_present(group_medians)),
                                 tar_target(lod_graph,
                                            graph_median_min(correlate_medians)))
  
  # running a single core to measure performance aspects ----------
performance_plan = tar_plan(
  single_core_perf = run_single_cor(),
  complexity_figure = create_complexity_figure(single_core_perf),
  
  egfrgenotype_n = c(1, 0.25, 0.5, 0.75, 0.99),
  tar_target(egfrgenotype_outliers,
             filter_generate_outliers(egfrgenotype_counts_info$counts, egfrgenotype_counts_info$info, egfrgenotype_n, "sample", "treatment"),
             pattern = map(egfrgenotype_n),
             iteration = "list"),
  egfrgenotype_single = get_single_outlier(egfrgenotype_outliers),

  tar_target(egfrgenotypetumorculture_outliers,
             filter_generate_outliers(egfrgenotypetumorculture_counts_info$counts, egfrgenotypetumorculture_counts_info$info, egfrgenotype_n, "sample", "treatment"),
             pattern = map(egfrgenotype_n),
             iteration = "list"),
  
  yeast_paper_outliers = c("WT.21", "WT.22", "WT.25", "WT.28", "WT.34", "WT.36",
                           "Snf2.06", "Snf2.13", "Snf2.25", "Snf2.35"),
  
  yeast_n = c(1, 0.25, 0.5),
  tar_target(yeast_outliers,
             filter_generate_outliers(yeast_counts_info$counts, yeast_counts_info$info, yeast_n, "sample", "treatment"),
             pattern = map(yeast_n),
             iteration = "list"),
  yeast_single = get_single_outlier(yeast_outliers),
  
  tar_target(adenocarcinoma_outliers,
             filter_generate_outliers(adenocarcinoma_counts_info$counts, adenocarcinoma_counts_info$info, yeast_n, "sample", "treatment"),
             pattern = map(yeast_n),
             iteration = "list"),
  adenocarcinoma_single = get_single_outlier(adenocarcinoma_outliers),
  
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
  
  tar_render(supp_materials,
              "doc/supplemental_materials.Rmd"),
  tar_render(supp_tables,
             "doc/supplemental_tables.Rmd"),
  tar_render(manuscript,
           "doc/ici_kt_manuscript.Rmd")

)

list(small_realistic_examples,
     loading_real_data,
     limit_of_detection_map,
     performance_plan)
