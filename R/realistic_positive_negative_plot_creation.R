compare_realistic_to_reference = function(realistic_sample_1,
                                      realistic_sample_2,
                                      realistic_neg_sample_2,
                                      realistic_na,
                                      realistic_positive_pearson,
                                      realistic_positive_kendall,
                                      realistic_positive_kt,
                                      realistic_negative_pearson_2,
                                      realistic_negative_kendall,
                                      realistic_negative_kt_2){
  

# tar_load(realistic_sample_1)
# tar_load(realistic_sample_2)
# tar_load(realistic_neg_sample_2)
  n_na = purrr::map_int(realistic_na, length)
  ref_pearson = cor(realistic_sample_1, realistic_sample_2)
ref_kendall = ici_kt(realistic_sample_1, realistic_sample_2, "global")[1]

ref_pearson_neg = cor(realistic_sample_1, realistic_neg_sample_2)
ref_kendall_neg = ici_kt(realistic_sample_1, realistic_neg_sample_2, "global")[1]

# tar_load(realistic_na)
# tar_load(realistic_positive_pearson)
# tar_load(realistic_positive_kendall)
# tar_load(realistic_positive_kt)
# tar_load(realistic_negative_pearson_2)
# tar_load(realistic_negative_kendall)
# tar_load(realistic_negative_kt_2)

realistic_positive_pearson = realistic_positive_pearson %>%
  dplyr::mutate(type = "Pearson",
                dir = "positive",
                diff = ref_pearson - cor,
                n_na = (x_na + y_na) / 2)

realistic_negative_pearson_2 = realistic_negative_pearson_2 %>%
  dplyr::mutate(type = "Pearson",
                dir = "negative",
                diff = ref_pearson_neg - cor,
                n_na = (x_na + y_na) / 2)

realistic_positive_kendall = realistic_positive_kendall %>%
  dplyr::mutate(type = "Kendall",
                dir = "positive",
                diff = ref_kendall - cor,
                n_na = (x_na + y_na) / 2)

realistic_negative_kendall = realistic_negative_kendall %>%
  dplyr::mutate(type = "Kendall",
                dir = "negative",
                diff = ref_kendall_neg - cor,
                n_na = (x_na + y_na) / 2)

realistic_positive_kt = realistic_positive_kt %>%
  dplyr::mutate(type = "ICI-Kt",
                dir = "positive",
                diff = ref_kendall - cor,
                n_na = (x_na + y_na) / 2)

realistic_negative_kt_2 = realistic_negative_kt_2 %>%
  dplyr::mutate(type = "ICI-Kt",
                dir = "negative",
                diff = ref_kendall_neg - cor,
                n_na = (x_na + y_na) / 2)

positive_df = rbind(realistic_positive_pearson,
                    realistic_positive_kendall,
                    realistic_positive_kt)

negative_df = rbind(realistic_negative_pearson_2,
                    realistic_negative_kendall,
                    realistic_negative_kt_2)

use_rows = rep(FALSE, nrow(negative_df))
n_subset = 10000
use_rows[sample(length(use_rows), n_subset)] = TRUE
negative_df = negative_df %>%
  dplyr::mutate(perc_missing = n_na / max(n_na) * 100,
                type = factor(type, levels = c("Kendall", "Pearson", "ICI-Kt"), ordered = TRUE)) %>%
  dplyr::filter(use_rows)
positive_df = positive_df %>%
  dplyr::mutate(perc_missing = n_na / max(n_na) * 100,
                type = factor(type, levels = c("Kendall", "Pearson", "ICI-Kt"), ordered = TRUE)) %>%
  dplyr::filter(use_rows)
pos_diff = ggplot(positive_df, aes(x = n_na, y = diff)) + 
  geom_point() +
  facet_grid(type ~ ., scales = "free") +
  labs(subtitle = "Positive Correlation",
       x = "Number Missing",
       y = "Difference from Reference") +
  scale_y_continuous(labels = label_scientific(digits = 2)) +
  theme(plot.subtitle = element_text(size = 13))
neg_diff = ggplot(negative_df, aes(x = n_na, y = diff)) +
  geom_point() +
  facet_grid(type ~ ., scales = "free") +
  labs(subtitle = "Negative Correlation",
       x = "Number Missing",
       y = "Difference from Reference") +
  scale_y_continuous(labels = label_scientific(digits = 2)) +
  theme(plot.subtitle = element_text(size = 13))

 list(positive = pos_diff,
      negative = neg_diff)
}

compare_realistic_to_each = function(realistic_sample_1,
                                          realistic_sample_2,
                                          realistic_neg_sample_2,
                                          realistic_na,
                                          realistic_positive_pearson,
                                          realistic_positive_kendall,
                                          realistic_positive_kt,
                                          realistic_negative_pearson_2,
                                          realistic_negative_kendall,
                                          realistic_negative_kt_2){
  
  
  # tar_load(realistic_sample_1)
  # tar_load(realistic_sample_2)
  # tar_load(realistic_neg_sample_2)
  # tar_load(realistic_na)
  # tar_load(realistic_positive_pearson)
  # tar_load(realistic_positive_kendall)
  # tar_load(realistic_positive_kt)
  # tar_load(realistic_negative_pearson_2)
  # tar_load(realistic_negative_kendall)
  # tar_load(realistic_negative_kt_2)
  
  # ref_pearson = cor(realistic_sample_1, realistic_sample_2)
  # ref_kendall = ici_kt(realistic_sample_1, realistic_sample_2, "global")[1]
  # 
  # ref_pearson_neg = cor(realistic_sample_1, realistic_neg_sample)
  # ref_kendall_neg = ici_kt(realistic_sample_1, realistic_neg_sample, "global")[1]
  
  
  n_na = purrr::map_int(realistic_na, length)
  
  rp_pearson_wide = realistic_positive_pearson |>
    dplyr::transmute(pearson = cor, id = paste0(i_na, ".", x_na, ".", y_na))
  
  rp_kendall_wide = realistic_positive_kendall |>
    dplyr::transmute(kendall = cor, id = paste0(i_na, ".", x_na, ".", y_na))
  
  rp_icikt_wide = realistic_positive_kt |>
    dplyr::transmute(icikt = cor, id = paste0(i_na, ".", x_na, ".", y_na))
  
  compare_positive = dplyr::left_join(rp_pearson_wide, rp_kendall_wide, by = "id")
  compare_positive = dplyr::left_join(compare_positive, rp_icikt_wide, by = "id")
  
  use_rows = rep(FALSE, nrow(compare_positive))
  n_subset = 10000
  use_rows[sample(length(use_rows), n_subset)] = TRUE
  
  compare_positive_subset = compare_positive[use_rows, ]
  
  rp_pearson_negative_wide = realistic_negative_pearson_2 %>%
    dplyr::transmute(pearson = cor,
                     id = paste0(i_na, ".", x_na, ".", y_na))
  rp_kendall_negative_wide = realistic_negative_kendall |>
    dplyr::transmute(kendall = cor,
                     id = paste0(i_na, ".", x_na, ".", y_na))
  rp_icikt_negative_wide = realistic_negative_kt_2 |>
    dplyr::transmute(icikt = cor,
                     id = paste0(i_na, ".", x_na, ".", y_na))
  
  compare_negative = dplyr::left_join(rp_pearson_negative_wide, rp_kendall_negative_wide, by = "id")
  compare_negative = dplyr::left_join(compare_negative, rp_icikt_negative_wide, by = "id")
  
  compare_negative_subset = compare_negative[use_rows, ]
  
  p_ici_pos = compare_positive_subset |>
    ggplot(aes(x = pearson, y = icikt)) + 
    geom_point() +
    labs(x = "Pearson", y = "ICI-Kt", subtitle = "Positive Correlation")
  k_ici_pos = compare_positive_subset |>
    ggplot(aes(x = kendall, y = icikt)) +
    geom_point() +
    labs(x = "Kendall", y = "ICI-Kt")
  
  p_ici_neg = compare_negative_subset |>
    ggplot(aes(x = pearson, y = icikt)) +
    geom_point() +
    labs(x = "Pearson", y = "ICI-Kt", subtitle = "Negative Correlation")
  k_ici_neg = compare_negative_subset |>
    ggplot(aes(x = kendall, y = icikt)) +
    geom_point() +
    labs(x = "Kendall", y = "ICI-Kt")
  
  
  
  list(positive = list(pearson = p_ici_pos,
                       kendall = k_ici_pos),
       negative = list(pearson = p_ici_neg,
                       kendall = k_ici_neg))
}
