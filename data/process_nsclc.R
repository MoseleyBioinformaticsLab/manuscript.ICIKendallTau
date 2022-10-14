# why?
# The NSCLC set of unassigned peaks across two instruments and 169 samples make
# a hefty dataset, 300 MB. We don't want it sitting in the _drake cache like
# that, so we are going to do the preprocessing and then load it in.
# 
library("ivs")
library(smirfeTools)
library(progress)
freq_to_iv = function(in_freq, in_sd, multiplier = 3)
{
  width_freq = in_sd * multiplier
  iv(start = in_freq - width_freq, end = in_freq + width_freq)
}

run_samples = function(in_sample)
{
  cluster1_peaks = all_peaks |>
    dplyr::filter(sample %in% in_sample)
  
  out_cluster1 = tibble::tibble(loc = vector("double", nrow(cluster1_peaks)),
                                sd = 0,
                                iv = iv(rep(1, nrow(cluster1_peaks)), rep(2, nrow(cluster1_peaks))),
                                peaklist = vector("list", nrow(cluster1_peaks)))
  
  save_index = 1
  pb = progress_bar$new(
    format = "   matching [:bar] :percent eta: :eta",
    total = nrow(cluster1_peaks),
    clear = FALSE)
  while (sum(!cluster1_peaks$used) > 0) {
    unused_peaks = cluster1_peaks |>
      dplyr::filter(!used)
    grab_peak1 = unused_peaks |>
      dplyr::slice_head(n = 1)
    matched_locations = iv_locate_overlaps(grab_peak1$iv, unused_peaks$iv)
    matched_peaks = unused_peaks[matched_locations$haystack, ]
    updated_peak1 = tibble::tibble(loc = median(matched_peaks$ObservedFrequency),
                                   sd = max(matched_peaks$ObservedFrequencySD),
                                   iv = freq_to_iv(loc, sd),
                                   peaklist = list(matched_peaks$PeakID))
    iteration_match = FALSE
    while (!iteration_match) {
      matched_locations2 = iv_locate_overlaps(updated_peak1$iv, unused_peaks$iv)
      iteration_match = length(base::setdiff(matched_locations2$haystack, matched_locations$haystack)) == 0
      matched_peaks = unused_peaks[matched_locations$haystack, ]
      updated_peak1 = tibble::tibble(loc = median(matched_peaks$ObservedFrequency),
                                     sd = max(matched_peaks$ObservedFrequencySD),
                                     iv = freq_to_iv(loc, sd),
                                     peaklist = list(matched_peaks$PeakID))
      matched_locations = matched_locations2
    }
    
    cluster1_peaks[cluster1_peaks$PeakID %in% updated_peak1$peaklist[[1]], "used"] = TRUE
    out_cluster1[save_index, ] = updated_peak1
    save_index = save_index + 1
    pb$update(ratio = sum(cluster1_peaks$used) / nrow(cluster1_peaks))
  }
  
  out_cluster1 = out_cluster1 |>
    dplyr::filter(!(sd %in% 0))
  out_cluster1
}



all_zip = dir("/big_data/data/nsclc_scpc_data/lung_matched_tissue-2022-05-11", pattern = "zip$", full.names = TRUE)
nsclc_data = readRDS("data/nsclc_peaks.rds")
names(nsclc_data) = purrr::map_chr(nsclc_data, \(.x) .x$Sample)

all_models = extract_coefficient_data(all_zip)
names(all_models) = gsub(".zip$", "", basename(all_zip))
all_coefficients = purrr::map_dfr(all_models, function(.x){
  tibble::tibble(coefficient = .x$coefficients$frequency_coefficients[3])
})
all_coefficients$sample = names(all_models)

all_coefficients = all_coefficients %>%
  dplyr::mutate(cluster =
                  dplyr::case_when(
                    coefficient < 29801700 ~ 1,
                    coefficient > 29801700 ~ 2
                  ))

all_peaks = purrr::imap(nsclc_data, function(.x, .y) {
  .x$Peaks |>
    dplyr::filter(!HighSD) |>
    dplyr::select(ObservedFrequency, ObservedFrequencySD, PeakID) |>
    dplyr::mutate(sample = .y,
                  PeakID = paste0(.y, ":", PeakID))
  })

all_peaks = tibble::as_tibble(dplyr::bind_rows(all_peaks))
all_peaks = dplyr::slice_sample(all_peaks, n = nrow(all_peaks))
all_peaks$used = FALSE
all_peaks$iv = freq_to_iv(all_peaks$ObservedFrequency, all_peaks$ObservedFrequencySD)

clustered_samples = split(all_coefficients$sample, all_coefficients$cluster)

clustered_peaks = purrr::map(clustered_samples, run_samples)



full_peak_list = tibble::tibble(mz = vector("double", n_peaks),
                                iv = iv(rep(1, n_peaks), rep(2, n_peaks)),
                                use = FALSE)


for (isample in all_peaks) {
  isample = isample |>
    dplyr::filter(!HighSD)
  isample_iv = tibble::tibble(mz = isample$ObservedMZ,
                              iv = mz_to_iv(isample$ObservedMZ),
                              use = TRUE)
  if (sum(full_peak_list$use) == 0) {
    full_peak_list[1:nrow(isample_iv), ] = isample_iv
  } else {
    tmp_full = full_peak_list |>
      dplyr::filter(use)
    aligned = iv_locate_overlaps(isample_iv$iv, tmp_full$iv)
    no_isample = aligned |>
      dplyr::filter(is.na(haystack))
    new_isample = isample_iv[no_isample$needles, ]
    start_peak = which.min(full_peak_list$use)
    end_peak = start_peak + nrow(new_isample) - 1
    message(sum(c(start_peak, length(start_peak:end_peak))))
    full_peak_list[start_peak:end_peak, ] = new_isample
  }
}

full_peak_list = full_peak_list[full_peak_list$use, ]

full_peak_list = full_peak_list |>
  dplyr::arrange(mz)

peak_matrix = matrix(NA, nrow = nrow(full_peak_list), ncol = length(all_peaks))
colnames(peak_matrix) = names(all_peaks)

for (isample in names(all_peaks))
{
  tmp_peaks = all_peaks[[isample]] |>
    dplyr::filter(!HighSD)
  tmp_iv = mz_to_iv(tmp_peaks$ObservedMZ)
  
  aligned = iv_locate_overlaps(full_peak_list$iv, tmp_iv) |>
    dplyr::filter(!is.na(haystack))
  peak_matrix[aligned$needles, isample] = tmp_peaks$Height[aligned$haystack]
}

rownames(peak_matrix) = paste0("f", seq(1, nrow(peak_matrix)))
saveRDS(peak_matrix, file = "data/nsclc_ppm_matched_peaks.rds")
