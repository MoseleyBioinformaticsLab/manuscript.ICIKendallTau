select_functions = rlang::syms(c("var_select", "pca_select"))
#reference_cor = readRDS(here::here("data/recount_adeno_cor.rds"))
fractions = (c(seq(0.01, 0.09, 0.01), seq(0.1, 0.9, 0.05)))
rep_10 = rep(0.1, 20)

small_samples = seq(4, 40, 4)
big_samples = seq(40, 264, 20)

set.seed(1234)

the_plan <-
  drake_plan(
   #set.seed(1234),
   x = seq(1, 10),
   y = seq(1, 10),
   y2 = seq(10, 1),
   where_na = create_na_indices(20),
   
   all_kt = all_kendalltau(x, y, y2, where_na),
   
   positive_kt = compare_positive_kt(x, y, where_na),
   negative_kt = compare_negative_kt(x, y2, where_na),
   positive_local_kt = compare_positive_kt(x, y, where_na, perspective = "local"),
   negative_local_kt = compare_negative_kt(x, y2, where_na, perspective = "local"),
   
   positive_low_kt = compare_positive_kt(x, y, where_na, low_indices = TRUE),
   negative_low_kt = compare_negative_kt(x, y2, where_na, low_indices = TRUE),
   positive_low_local_kt = compare_positive_kt(x, y, where_na, low_indices = TRUE, perspective = "local"),
   negative_low_local_kt = compare_negative_kt(x, y2, where_na, low_indices = TRUE, perspective = "local"),
   
   positive_pearson = compare_positive_pearson(x, y, where_na),
   negative_pearson = compare_negative_pearson(x, y2, where_na),
   
   positive_low_pearson = compare_positive_pearson(x, y, where_na, low_indices = TRUE),
   negative_low_pearson = compare_negative_pearson(x, y2, where_na, low_indices = TRUE),
   
   positive_low_kendall = compare_positive_pearson(x, y, where_na, low_indices = TRUE, method = "kendall"),
   negative_low_kendall = compare_negative_pearson(x, y2, where_na, low_indices = TRUE, method = "kendall"),
   
   positive_kendall = compare_positive_pearson(x, y, where_na, method = "kendall"),
   negative_kendall = compare_negative_pearson(x, y2, where_na, method = "kendall"),
   
   realistic_sample_1 = create_sample(n = 1000),
   realistic_sample_2 = create_sample(n = 1000),
   realistic_neg_sample = sort(realistic_sample_2, decreasing = TRUE),
   realistic_na = create_random_na(),
   realistic_positive_kt = compare_positive_kt(realistic_sample_1, realistic_sample_2, realistic_na),
   realistic_negative_kt = compare_negative_kt(realistic_sample_1, realistic_neg_sample, realistic_na),
   
   realistic_positive_pearson = compare_positive_pearson(realistic_sample_1, realistic_sample_2, realistic_na),
   realistic_negative_pearson = compare_negative_pearson(realistic_sample_1, realistic_neg_sample, realistic_na),
   
   realistic_positive_kendall = compare_positive_pearson(realistic_sample_1, realistic_sample_2, realistic_na, method = "kendall"),
   realistic_negative_kendall = compare_negative_pearson(realistic_sample_1, realistic_neg_sample, realistic_na, method = "kendall"),
   
   # kendall_pearson_comparison = target(
   #   command = {
   #     rmarkdown::render(knitr_in("doc/kendall_pearson_comparisons.Rmd"))
   #     file_out("doc/kendall_pearson_comparisons.pdf")
   #   }
   # ),
   
   transcript_data = readRDS(here::here("data/recount_adeno_counts.rds")),
   transcript_pca = prcomp(t(log1p(transcript_data)), center = TRUE, scale. = FALSE),
   
   transcript_na = zero_to_na(transcript_data),
   
   ref_cor = ici_kendalltau(t(transcript_na)),
   ref_completeness = pairwise_completeness(t(transcript_na)),
   
   select_nonrandom_fraction = target(
      select_function(transcript_pca, transcript_na, fraction = frac_value),
      transform = cross(
         frac_value = !!fractions,
         select_function = !!select_functions
      )
   ),
   
   run_nonrandom = target(
      run_fractional_correlation(select_nonrandom_fraction),
      transform = map(select_nonrandom_fraction)
   ),
   
   select_random_fraction = target(
      select_random(transcript_na, fraction = frac_value),
      transform = map(frac_value = !!fractions)
   ),
   
   run_random = target(
      run_fractional_correlation(select_random_fraction),
      transform = map(select_random_fraction)
   ),
   
   # compare_nonrandom = target(
   #    compare_fractional_correlation(run_nonrandom_fraction, reference_cor),
   #    transform = map(run_nonrandom_fraction)
   # ),
   # 
   # compare_random = target(
   #    compare_fractional_correlation(run_random_fraction, reference_cor),
   #    transform = map(run_random_fraction)
   # ),
   
   select_random_multiple = target(
      select_random(transcript_na, fraction = frac_value),
      transform = map(frac_value = !!rep_10)
   ),
   
   run_random_multiple = target(
      run_fractional_correlation(select_random_multiple),
      transform = map(select_random_multiple)
   ),
   
   select_ss_small = target(
      select_samples(transcript_na, n_sample = n_sample),
      transform = map(n_sample = !!small_samples)
   ),
   
   select_ss_big = target(
      select_samples(transcript_na, n_sample = n_sample),
      transform = map(n_sample = !!big_samples)
   ),
   
   run_small = target(
      run_small_samples(select_ss_small),
      transform = map(select_ss_small)
   ),
   
   run_big = target(
     run_big_samples(select_ss_big),
     transform = map(select_ss_big)
   ),
   
   results_random = target(
      random_2_reference(run_random, ref_cor),
      transform = map(run_random)
   ),
   
   combined_random = target(
      bind_rows(results_random),
      transform = combine(results_random)
   ),
   
   results_small = target(
      get_run_time(run_small),
      transform = map(run_small)
   ),
   
   combined_small = target(
      bind_rows(results_small),
      transform = combine(results_small)
   ),
   
   results_big = target(
      get_run_time(run_big),
      transform = map(run_big)
   ),
   
   combined_big = target(
      bind_rows(results_big),
      transform = combine(results_big)
   ),
   
   results_nonrandom = target(
      random_2_reference(run_nonrandom, ref_cor),
      transform = map(run_nonrandom)
   ),
   
   combined_nonrandom = target(
      bind_rows(results_nonrandom),
      transform = combine(results_nonrandom)
   ),
   
   results_rand_multiple = target(
      random_2_reference(run_random_multiple, ref_cor),
      transform = map(run_random_multiple)
   ),
   
   combined_rand_multiple = target(
      bind_rows(results_rand_multiple),
      transform = combine(results_rand_multiple)
   ),
   
   reshaped_random = reshape_data(combined_random),
   reshaped_nonrandom = reshape_data(combined_nonrandom),
   reshaped_multiple = reshape_data(combined_rand_multiple),
   
   # pearson and kendall results
   setup_pearson = list(data = transcript_na,
                        type = "random",
                        frac = 1),
   ref_pearson = run_fractional_pearson(setup_pearson),
   ref_pearson_0 = run_fractional_pearson(setup_pearson, replace_0 = TRUE),
   ref_kendall = run_fractional_kendall(setup_pearson),
   ref_kendall_0 = run_fractional_kendall(setup_pearson, replace_0 = TRUE),
   
   run_random_pearson = target(
      run_fractional_pearson(select_random_fraction),
      transform = map(select_random_fraction)
   ),
   
   run_random_pearson_0 = target(
      run_fractional_pearson(select_random_fraction, replace_0 = TRUE),
      transform = map(select_random_fraction)
   ),
   
   run_random_kendall = target(
      run_fractional_kendall(select_random_fraction),
      transform = map(select_random_fraction)
   ),
   
   run_random_kendall_0 = target(
      run_fractional_kendall(select_random_fraction, replace_0 = TRUE),
      transform = map(select_random_fraction)
   ),
   
   left_censored_samples = create_lc_samples(),
   left_censored_cor = left_censor_correlate(left_censored_samples),
   random_censored_cor = random_censor_correlate(left_censored_samples),
   logtransform_censored_cor = lt_left_censor_correlate(left_censored_samples),
   
   
   brainsonrnaseq_counts = readRDS(here::here("data", "brainson_rnaseq201901_counts.rds")),
   brainsonrnaseq_info = readRDS(here::here("data", "brainson_rnaseq201901_info.rds")),
   brainsonrnaseq_outliers_1 = filter_generate_outliers(brainsonrnaseq_counts, brainsonrnaseq_info, 1, "sample", "tumor"),
   brainsonrnaseq_outliers_25 = filter_generate_outliers(brainsonrnaseq_counts, brainsonrnaseq_info, 0.25, "sample", "tumor"),
   brainsonrnaseq_outliers_50 = filter_generate_outliers(brainsonrnaseq_counts, brainsonrnaseq_info, 0.5, "sample", "tumor"),
   
   yeast_paper_outliers = c("WT.21", "WT.22", "WT.25", "WT.28", "WT.34", "WT.36",
                            "SNF2.06", "SNF2.13", "SNF2.25", "SNF2.35"),
   yeast_counts_info = readRDS(here::here("data", "yeast_counts_info.rds")),
   yeast_outliers_1 = filter_generate_outliers(yeast_counts_info$counts, yeast_counts_info$info, 1,
                                               "sample_rep", "sample"),
   yeast_outliers_50 = filter_generate_outliers(yeast_counts_info$counts, yeast_counts_info$info, 0.5,
                                               "sample_rep", "sample"),
   yeast_outliers_25 = filter_generate_outliers(yeast_counts_info$counts, yeast_counts_info$info, 0.25,
                                               "sample_rep", "sample"),
   
   adeno_info = readRDS(here::here("data", "transcript_info.rds")),
   adeno_data = transcript_data,
   adeno_outliers_1 = filter_generate_outliers(adeno_data, adeno_info, 1,
                                               "sample_id2", "tissue_type"),
   adeno_outliers_25 = filter_generate_outliers(adeno_data, adeno_info, 0.25,
                                               "sample_id2", "tissue_type"),
   adeno_outliers_50 = filter_generate_outliers(adeno_data, adeno_info, 0.5,
                                               "sample_id2", "tissue_type"),
   
   eval_random = target(
      evaluate_by_pca(select_random_fraction, transcript_pca),
      transform = map(select_random_fraction)
   ),
   
   combined_eval_random = target(
      bind_rows(eval_random),
      transform = combine(eval_random)
   ),
   
   eval_nonrandom = target(
      evaluate_by_pca(select_nonrandom_fraction, transcript_pca),
      transform = map(select_nonrandom_fraction)
   ),
   
   combined_eval_nonrandom = target(
      bind_rows(eval_nonrandom),
      transform = combine(eval_nonrandom)
   ),
   
   all_pca_eval = bind_rows(combined_eval_random,
                            combined_eval_nonrandom),
   
   pca_eval_summary = summarize_by_pca(all_pca_eval),
   
   single_core_perf = run_single_cor(),
   complexity_figure = create_complexity_figure(single_core_perf)

)
