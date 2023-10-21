library(ivs)
imf_data = readRDS(here::here("raw_data/nsclc/nsclc_scancentric_imfs_raw"))
class_data = smirfeTools::import_emf_classifications(file.path(here::here("raw_data/nsclc"), "all_emfs_classified_2020-01-27.json"))

other_data = readRDS(here::here("data/nsclc_other_peak_data.rds"))

mz_to_iv = function(in_mz, in_ppm)
{
  width_mz = in_ppm * 1e-6 * in_mz
  iv(start = in_mz - width_mz, end = in_mz + width_mz)
}


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

all_samples = unique(other_data$all_peaks$sample)
base::setdiff(all_samples, colnames(imf_peak_lists))


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

matched_ids = dplyr::left_join(matched_ids, 
