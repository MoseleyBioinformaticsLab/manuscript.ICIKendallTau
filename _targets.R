## Load your packages, e.g. library(targets). -----
tar_source("./packages.R")

## Load your R files -------
tar_source(files = "R")

# creating mapping tibble --------
dataset_variables = tibble::tibble(character = c("yeast_counts_info",
                                                 "adenocarcinoma_counts_info",
                                  "egfrgenotype_counts_info",
                                  "egfrgenotypetumorculture_counts_info",
                                  "typeandtumorculture_counts_info",
                                  "ratstamina_counts_info",
                                  "nsclc_counts_info")) |>
  dplyr::mutate(id = gsub("_.*", "", character),
                sym = rlang::syms(character),
                Measurement = c("RNA-Seq",
                                "RNA-Seq",
                                "RNA-Seq",
                                "RNA-Seq",
                                "RNA-Seq",
                                "Metabolomics",
                                "Lipidomics"))

correlation_methods = c("ici",
                                    "ici_completeness",
                                    "pearson_base",
                                    "pearson_base_nozero",
                                    "pearson_log1p",
                                    "pearson_log",
                                    "kt")

dataset_feature_correlation = tidyr::expand_grid(character_dataset = dataset_variables$sym,
                                                 character_correlation = correlation_methods) |>
  dplyr::mutate(dataset = rlang::syms(character_dataset),
                correlation = rlang::syms(character_correlation),
                 id = paste0(character_correlation, "_", gsub("_.*", "", character_dataset)))
  # dplyr::filter(grepl("yeast|ratstamina", id))

# realistic examples -----
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
  realistic_neg_sample_2 = -1 * realistic_sample_2,
  realistic_na = create_random_na(),
  realistic_positive_kt = compare_positive_kt(realistic_sample_1, realistic_sample_2, realistic_na),
  realistic_negative_kt = compare_negative_kt(realistic_sample_1, realistic_neg_sample, realistic_na),
  realistic_negative_kt_2 = compare_negative_kt(realistic_sample_1, realistic_neg_sample_2, realistic_na),
  
  realistic_positive_pearson = compare_positive_pearson(realistic_sample_1, realistic_sample_2, realistic_na),
  realistic_negative_pearson = compare_negative_pearson(realistic_sample_1, realistic_neg_sample, realistic_na),
  realistic_negative_pearson_2 = compare_negative_pearson(realistic_sample_1, realistic_neg_sample_2, realistic_na),
  
  realistic_positive_kendall = compare_positive_kendall(realistic_sample_1, realistic_sample_2, realistic_na),
  realistic_negative_kendall = compare_negative_kendall(realistic_sample_1, realistic_neg_sample, realistic_na),
  
  realistic_reference_plot = compare_realistic_to_reference(realistic_sample_1,
                                                         realistic_sample_2,
                                                         realistic_neg_sample_2,
                                                         realistic_na,
                                                         realistic_positive_pearson,
                                                         realistic_positive_kendall,
                                                         realistic_positive_kt,
                                                         realistic_negative_pearson_2,
                                                         realistic_negative_kendall,
                                                         realistic_negative_kt_2),
  realistic_comparison_plot = compare_realistic_to_each(realistic_sample_1,
                                                        realistic_sample_2,
                                                        realistic_neg_sample_2,
                                                        realistic_na,
                                                        realistic_positive_pearson,
                                                        realistic_positive_kendall,
                                                        realistic_positive_kt,
                                                        realistic_negative_pearson_2,
                                                        realistic_negative_kendall,
                                                        realistic_negative_kt_2),
  
  left_censored_samples = create_lc_samples(),
  left_censored_cor = left_censor_correlate(left_censored_samples),
  random_censored_cor = random_censor_correlate(left_censored_samples),
  logtransform_censored_cor = lt_left_censor_correlate(left_censored_samples),

  subsample = sample(1000, 50),
  left_sample_cor = left_censor_correlate(left_censored_samples[subsample, ]),
  random_sample_cor = random_censor_correlate(left_censored_samples[subsample, ], n_na = seq(0, 12, 2)),
  logtransform_sample_cor = lt_left_censor_correlate(left_censored_samples[subsample, ]),
  
  censored_value_plots = plot_censored_data(left_censored_cor,
                                            random_censored_cor,
                                            logtransform_censored_cor),
  censored_compare_plots = compare_censored_data(left_censored_cor,
                                                 random_censored_cor,
                                                 logtransform_censored_cor)
)



# loading real data ---------
loading_real_data = tar_plan(
  tar_target(adenocarcinoma_file,
             "data/recount_adenocarcinoma_count_info.rds",
             format = "file"),
  adenocarcinoma_counts_info = readRDS(adenocarcinoma_file),
  
  tar_target(barton_yeast_file,
             "data/barton_yeast_counts_info.rds",
             format = "file"),
  yeast_counts_info = readRDS(barton_yeast_file),
  
  tar_target(brainson_egfrgenotype_file,
             "data/brainsonrnaseq_type_counts_info.rds",
             format = "file"),
  egfrgenotype_counts_info = readRDS(brainson_egfrgenotype_file),
  
  tar_target(brainson_tumorcultures_file,
             "data/brainsonrnaseq_type_counts_info_tumor.rds",
             format = "file"),
  egfrgenotypetumorculture_counts_info = readRDS(brainson_tumorcultures_file),
  
  tar_target(brainson_typeandtumorcultures_file,
             "data/brainsonrnaseq_type_counts_info_typetumor.rds",
             format = "file"),
  typeandtumorculture_counts_info = readRDS(brainson_typeandtumorcultures_file),
  
  tar_target(mwtab_ratstamina_file,
             "data/mwtab_st000017_an000034_count_info.rds",
             format = "file"),
  ratstamina_counts_info = readRDS(mwtab_ratstamina_file),
  
  tar_target(rcsirm_nsclc_file,
             "data/nsclc_count_info.rds",
             format = "file"),
  nsclc_counts_info = readRDS(rcsirm_nsclc_file)
)

# limit of detection examples --------


limit_of_detection_map = tar_map(dataset_variables,
                                    names = id,
                                    tar_target(rank_ordered,
                                               rank_order_datasets(sym, id)),
                                    tar_target(correlate_ranks,
                                      calculate_median_min_correlation(rank_ordered)),
                                    tar_target(rank_graphs,
                                               graph_median_min(correlate_ranks)))

# left-censorship calculation ------
left_censorship_map = tar_map(dataset_variables,
                              names = id,
                              tar_target(lc_test,
                                         test_left_censorship_datasets(sym, id)),
                              tar_target(lc_summary,
                                         create_lc_summary(lc_test)))


left_censorship_combine = tar_combine(lc_summary_all,
                                      left_censorship_map[[2]],
                                      command = bind_rows(!!!.x))
  
# running a single core to measure performance aspects ----------
performance_plan = tar_plan(
  single_core_perf = run_single_cor(),
  complexity_figure = create_complexity_figure(single_core_perf),
  
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
             pattern = map(run_multi))
)

# sample outliers -----
sample_outlier_plan = tar_plan(
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
  
  egfrgenotypetumorculture_single = get_single_outlier(egfrgenotypetumorculture_outliers),
  
  tar_target(typeandtumorculture_outliers,
             filter_generate_outliers(typeandtumorculture_counts_info$counts,
                                      typeandtumorculture_counts_info$info,
                                      egfrgenotype_n, "sample", "treatment"),
             pattern = map(egfrgenotype_n),
             iteration = "list"),
  typeandtumorculture_single = get_single_outlier(typeandtumorculture_outliers),
  
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
  
  tar_target(nsclc_outliers,
             filter_generate_outliers(nsclc_counts_info$counts, nsclc_counts_info$info, yeast_n, "sample", "treatment"),
             pattern = map(yeast_n),
             iteration = "list"),
  nsclc_single = get_single_outlier(nsclc_outliers),
  
  tar_target(ratstamina_outliers,
             filter_generate_outliers(ratstamina_counts_info$counts, ratstamina_counts_info$info, yeast_n, "sample", "treatment"),
             pattern = map(yeast_n),
             iteration = "list"),
  ratstamina_single = get_single_outlier(ratstamina_outliers),
  
  ## Annotations -----
  tar_target(kegg_data,
             "data/kegg_compound_mapping.rds",
             format = "file"),
  tar_target(nsclc_lipidclasses,
             "data/nsclc_feature_lipid_classes.rds",
             format = "file"),
  feature_annotations = get_feature_annotations(kegg_data, nsclc_lipidclasses),
  tar_target(metabolite_kegg_file,
             "data/mwtab_st000017_an000034_metabolites.rds",
             format = "file"),
  metabolite_kegg = readRDS(metabolite_kegg_file)
  
)


# correlation of features -----
feature_correlation_map = tar_map(dataset_feature_correlation,
                                 names = id,
                                 tar_target(feature_correlation,
                                            correlation(dataset, id, 0.25, "sample", "treatment")),
                                 tar_target(feature_partial_cor,
                                            calculate_partial_cor_pvalues(feature_correlation)),
                                 tar_target(feature_network_qratio,
                                            calculate_feature_network_qratio(feature_partial_cor,
                                                                             feature_annotations, 
                                                                             "pathway", 
                                                                             metabolite_kegg)),
                                 tar_target(feature_network_qratio_w1,
                                            calculate_feature_network_qratio(feature_partial_cor,
                                                                             feature_annotations, 
                                                                             "pathway", 
                                                                             metabolite_kegg,
                                                                             use_weights = FALSE)))


feature_qratio_combine_map = tar_combine(feature_qratio_comparisons,
                                         feature_correlation_map[[3]],
                                         command = bind_rows(!!!.x))

feature_qratio_summary_plan = tar_plan(
  feature_qratio_summary = feature_qratio_comparisons |>
    dplyr::filter(!is.na(q_value)) |>
    dplyr::select(data_id, method_id, q_value) |>
    dplyr::distinct() |>
    dplyr::group_by(data_id) |>
    dplyr::arrange(dplyr::desc(q_value), .by_group = TRUE)
)




## dataset summaries -----
dataset_summary_map = tar_map(dataset_variables,
                               names = id,
                               tar_target(summary_n,
                                          create_dataset_summary(sym)))
dataset_summary_combine = tar_combine(dataset_summary_n,
                                      dataset_summary_map,
                                      command = bind_rows(!!!.x))

dataset_summary_plan = tar_plan(dataset_summary = dplyr::bind_cols(
  dataset_variables[, c("id", "Measurement")], dataset_summary_n
) |>
  dplyr::mutate(Dataset = id,
                id = NULL) |>
  dplyr::select(Dataset, Measurement, Features, Samples, Conditions, Replicates)
)

# documents -----
documents_plan = tar_plan(
  tar_target(references,
             "docs/icikt_references.json",
             format = "file"),
  tar_render(supp_materials,
              "docs/supplemental_materials.Rmd"),
  tar_render(manuscript,
           "docs/ici_kt_manuscript.Rmd"),

)

# put all together -----
list(small_realistic_examples,
     loading_real_data,
     limit_of_detection_map,
     left_censorship_map,
     left_censorship_combine,
     sample_outlier_plan,
     feature_correlation_map,
     feature_qratio_combine_map,
     feature_qratio_summary_plan,
     dataset_summary_map,
     dataset_summary_combine,
     dataset_summary_plan,
     performance_plan,
     documents_plan)
