## useful functions -----
mz_to_iv = function(in_mz, in_ppm)
{
  width_mz = in_ppm * 1e-6 * in_mz
  iv(start = in_mz - width_mz, end = in_mz + width_mz)
}

voted_2_multi = function(voted_categories){
  split_by_emf = split(voted_categories, voted_categories$emf)
  
  purrr::map_df(split_by_emf, function(.x){
    if (nrow(.x) == 1) {
      return(.x)
    } else {
      tmp_df = .x[1, ]
      tmp_df$Categories = "multiple"
      tmp_df$counts = nrow(.x)
      tmp_df$percent = 100
      return(tmp_df)
    }
  })
}

merge_voted_classes = function(imf_voted_categories, imf_peak_categories){
  imf_peak_categories_voted_classes = purrr::map_dfr(seq(1, nrow(imf_voted_categories)), function(in_row){
    tmp_data = imf_voted_categories[in_row, ]
    match_imf = dplyr::filter(imf_peak_categories, sudo_EMF %in% tmp_data$sudo_EMF)
    if (tmp_data$Categories %in% "multiple") {
      match_imf2 = dplyr::select(match_imf, -Categories) %>% unique()
    } else {
      match_imf2 = dplyr::filter(match_imf, Categories %in% tmp_data$Categories) %>%
        dplyr::select(-Categories) %>% unique()
    }
    match_imf2$VotedCategories = tmp_data$Categories
    match_imf2
  })
  imf_peak_categories_voted_classes
}

## data -----
library(ivs)
imf_data = readRDS(here::here("raw_data/nsclc/nsclc_scancentric_imfs_raw"))
imf_info = imf_data$info
class_data = smirfeTools::import_emf_classifications(file.path(here::here("raw_data/nsclc"), "all_emfs_classified_2020-01-27.json"))
emf_info = readRDS(here::here("raw_data/nsclc/nsclc_emf_info.rds"))
other_data = readRDS(here::here("data/nsclc_other_peak_data.rds"))


## mapping -----
# now we can grab just the complete to isotopologue EMF mapping that is necessary
# to go from the class data (isotopologue) to peak data (complete)
complete_isotopologue_mapping = emf_info |>
  dplyr::select(complete_EMF, isotopologue_EMF) |>
  dplyr::distinct()

# and map the classes to the peaks.
class_2_complete_emf = dplyr::left_join(class_data, complete_isotopologue_mapping,
                              by = "isotopologue_EMF",
                              relationship = "many-to-many")

## match peaks between data sets ----
imf_peak_lists = imf_data$peaks
imf_peak_lists = sub("_", ":", imf_peak_lists)
imf_peaks = purrr::map(seq_len(nrow(imf_peak_lists)), \(in_row){
  tmp = imf_peak_lists[in_row, ]
  names(tmp) = NULL
  tmp = tmp[!is.na(tmp)]
  tmp
})
matched_peaks = other_data$matched_peaks
matched_peaks$feature_id = paste0("f", seq_len(nrow(matched_peaks)))
matched_peaks$iv2 = iv(start = matched_peaks$loc, end = matched_peaks$loc + 1e-8)

imf_peak_locs = imf_data$location
imf_peak_df = tibble::tibble(peak = rownames(imf_peak_locs), location = rowMedians(imf_peak_locs, na.rm = TRUE))
imf_peak_df$peak_list = imf_peaks
imf_peak_df$iv = mz_to_iv(imf_peak_df$location, 0.5)

matches = iv_locate_overlaps(imf_peak_df$iv, matched_peaks$iv2)
matches = matches |>
  dplyr::filter(!is.na(haystack))

n_match = matches |>
  dplyr::group_by(needles) |>
  dplyr::summarise(n = dplyr::n())

single_match = n_match |>
  dplyr::filter(n == 1)

matches = matches |>
  dplyr::filter(needles %in% single_match$needles)

matched_ids = tibble::tibble(sudo_EMF = imf_peak_df$peak[matches$needles],
                             p_sudo = imf_peak_df$peak_list[matches$needles], 
                             feature_id = matched_peaks$feature_id[matches$haystack],
                             p_feature = matched_peaks$peaklist[matches$haystack])

n_peak_match = purrr::map_dbl(seq_len(nrow(matched_ids)), \(in_row){
  perc_match = length(base::intersect(matched_ids$p_sudo[in_row][[1]], matched_ids$p_feature[in_row][[1]])) / length(matched_ids$p_sudo[in_row][[1]])
  perc_match
})

matched_ids$perc_match = n_peak_match

matched_ids_trust = matched_ids |>
  dplyr::filter(perc_match >= 0.9)

matched_ids_trust = dplyr::left_join(matched_ids_trust, imf_data$info, by = c("sudo_EMF" = "PeakID"))

## add classes to feature_id ----
feature_id_2_class = dplyr::left_join(matched_ids_trust[, c("feature_id", "complete_EMF", "sudo_EMF")],
                                      class_2_complete_emf, by = c("complete_EMF"),
                                      relationship = "many-to-many")

## vote on them to make a decision! ----
trimmed_sudo_emf_class = feature_id_2_class |>
  dplyr::mutate(emf = gsub(".IMF_.*", "", sudo_EMF)) |>
  dplyr::select(emf, complete_EMF, Categories, Classes) |>
  dplyr::distinct()
voted_data = metabolomicsUtilities::vote_categories_classes(trimmed_sudo_emf_class)

feature_id_2_class = feature_id_2_class |>
    dplyr::mutate(emf = gsub(".IMF_.*", "", sudo_EMF))

feature_id_2_voted = dplyr::left_join(feature_id_2_class[, c("feature_id", "sudo_EMF", "emf")], voted_data,
                                      by = "emf")
saveRDS(feature_id_2_voted, file = "data/nsclc_feature_lipid_classes.rds")

