# why?
# The NSCLC set of unassigned peaks across two instruments and 169 samples make
# a hefty dataset, 300 MB. We don't want it sitting in the _drake cache like
# that, so we are going to do the preprocessing and then load it in.
# 
library(ivs)
library(progress)
library(ggplot2)
library(furrr)
options(parallelly.fork.enable = TRUE)
plan(multicore)
mz_to_iv = function(in_mz, in_ppm)
{
  width_mz = in_ppm * 1e-6 * in_mz
  iv(start = in_mz - width_mz, end = in_mz + width_mz)
}

match_peaks = function(in_peaks)
{
  save_index = 1
  out_peaks = tibble::tibble(loc = vector("double", nrow(in_peaks)),
                                sd = 0,
                                iv = iv(rep(1, nrow(in_peaks)), rep(2, nrow(in_peaks))),
                                peaklist = vector("list", nrow(in_peaks)))
  while (sum(!in_peaks$used) > 0) {
    unused_peaks = in_peaks |>
      dplyr::filter(!used)
    grab_peak1 = unused_peaks |>
      dplyr::slice_head(n = 1)
    matched_locations = iv_locate_overlaps(grab_peak1$iv, unused_peaks$iv)
    matched_peaks = unused_peaks[matched_locations$haystack, ]
    updated_peak1 = tibble::tibble(loc = median(matched_peaks$ObservedMZ),
                                   sd = sd(matched_peaks$ObservedMZ),
                                   iv = mz_to_iv(loc, 0.5),
                                   peaklist = list(matched_peaks$PeakID))
    iteration_match = FALSE
    while (!iteration_match) {
      matched_locations2 = iv_locate_overlaps(updated_peak1$iv, unused_peaks$iv)
      iteration_match = length(base::setdiff(matched_locations2$haystack, matched_locations$haystack)) == 0
      matched_peaks = unused_peaks[matched_locations$haystack, ]
      updated_peak1 = tibble::tibble(loc = median(matched_peaks$ObservedMZ),
                                     sd = sd(matched_peaks$ObservedMZ),
                                     iv = mz_to_iv(loc, 0.5),
                                     peaklist = list(matched_peaks$PeakID))
      matched_locations = matched_locations2
    }
    
    in_peaks[in_peaks$PeakID %in% updated_peak1$peaklist[[1]], "used"] = TRUE
    out_peaks[save_index, ] = updated_peak1
    save_index = save_index + 1
  }
  
  out_peaks = out_peaks |>
    dplyr::filter(!(sd %in% 0))
  out_peaks
}

run_samples = function(run_peaks)
{
  run_peaks$mziv = iv(start = run_peaks$ObservedMZ, end = run_peaks$ObservedMZ + 1e-8)
  
  ms_ranges = iv(start = seq(150, 1600 - 50, by = 50), end = seq(200, 1600, by = 50))
  
  ms_overlaps = iv_locate_overlaps(ms_ranges, run_peaks$mziv)
  
  split_overlaps = split(ms_overlaps$haystack, ms_overlaps$needles)
  split_ms = purrr::map(split_overlaps, ~ run_peaks[.x, ])
  
  run_each = furrr::future_map(split_ms, match_peaks)
  
  matched_peaks = dplyr::bind_rows(run_each)
  matched_peaks
}

## work with assigned data first ----
## This is going to help us figure out the limit to use to match the peaks
## to each other.
imf_data = readRDS(here::here("raw_data/nsclc/nsclc_scancentric_imfs_raw"))
imf_loc_data = tibble::tibble(sd = apply(imf_data$location, 1, sd, na.rm = TRUE),
                              mean = apply(imf_data$location, 1, mean, na.rm = TRUE)) |>
  dplyr::mutate(sd_ppm = sd / mean / 1e-6) |>
  dplyr::filter(!is.na(sd), sd_ppm < 10)

imf_loc_data |>
  ggplot(aes(x = sd_ppm)) +
  geom_histogram(bins = 100)

## Based on this graph, a limit of 0.5 to either side seems like we would be OK.
nsclc_data = readRDS(here::here("raw_data/nsclc/nsclc_peaks.rds"))
names(nsclc_data) = purrr::map_chr(nsclc_data, \(.x) .x$Sample)

all_peaks = purrr::imap(nsclc_data, function(.x, .y) {
  .x$Peaks |>
    dplyr::filter(!HighSD) |>
    dplyr::select(ObservedMZ, ObservedMZSD, PeakID, Height) |>
    dplyr::mutate(sample = .y,
                  PeakID = paste0(.y, ":", PeakID))
  })

all_peaks = tibble::as_tibble(dplyr::bind_rows(all_peaks))
all_peaks$used = FALSE
all_peaks$iv = mz_to_iv(all_peaks$ObservedMZ, 0.5)

nsclc_info = readRDS(here::here("raw_data/nsclc/nsclc_info"))
match_samples = base::intersect(all_peaks$sample, nsclc_info$sample)

nsclc_info = nsclc_info |>
  dplyr::filter(sample %in% match_samples)
rownames(nsclc_info) = nsclc_info$sample

all_peaks = all_peaks |>
  dplyr::filter(sample %in% match_samples)
# checking that the number of peaks varies only a little
# if we scramble the data differently
# use_seeds = c(1234, 9034, 3042, 5467)
# 
# for (iseed in use_seeds) {
#   set.seed(iseed)
#   all_peaks = dplyr::slice_sample(all_peaks, n = nrow(all_peaks))
#   matched_peaks = run_samples(all_peaks)
#   message(nrow(matched_peaks))
# }

set.seed(5467)
all_peaks = dplyr::slice_sample(all_peaks, n = nrow(all_peaks))
matched_peaks = run_samples(all_peaks)

matched_heights = matrix(NA, nrow = nrow(matched_peaks), ncol = length(unique(all_peaks$sample)))
matched_locations = matrix(NA, nrow = nrow(matched_peaks), ncol = length(unique(all_peaks$sample)))
colnames(matched_heights) = colnames(matched_locations) = unique(all_peaks$sample)
pb = progress_bar$new(total = nrow(matched_heights),
                      format = ":spin [:bar] :percent in :elapsed ETA: :eta")
for (ipeak in seq_len(nrow(matched_heights))) {
  tmp_peaklist = matched_peaks$peaklist[ipeak][[1]]
  tmp_peaks = all_peaks |>
    dplyr::filter(PeakID %in% tmp_peaklist)
  matched_heights[ipeak, tmp_peaks$sample] = tmp_peaks$Height
  matched_locations[ipeak, tmp_peaks$sample] = tmp_peaks$ObservedMZ
  if ((ipeak %% 10) == 0) {
    pb$tick(10)
  }
}

class_data = smirfeTools::import_emf_classifications(file.path(here::here("raw_data/nsclc"), "all_emfs_classified_2020-01-27.json"))


nsclc_medians = readRDS(here::here("raw_data/nsclc/nsclc_scancentric_medians"))
nsclc_medians = nsclc_medians |>
  dplyr::filter(sample %in% match_samples)
median_matrix = matrix(nsclc_medians$median, nrow = nrow(matched_heights),
                       ncol = ncol(matched_heights), byrow = TRUE)
colnames(median_matrix) = nsclc_medians$sample
median_matrix = median_matrix[, colnames(matched_heights)]

matched_heights_norm = matched_heights / median_matrix
matched_heights_norm[is.na(matched_heights_norm)] = 0
rownames(matched_heights_norm) = paste0("f", seq_len(nrow(matched_heights_norm)))
nsclc_info = nsclc_info |>
  dplyr::mutate(treatment = paste0(instrument, ":", disease))

matched_heights_norm = matched_heights_norm[, match_samples]
nsclc_info = nsclc_info[match_samples, ]
all.equal(colnames(matched_heights_norm), nsclc_info$sample)


out_norm = list(counts = matched_heights_norm, info = nsclc_info)
out_norm$data_id = "nsclc"
saveRDS(out_norm, here::here("data", "nsclc_count_info.rds"))
other_data = list(matched_peaks = matched_peaks, all_peaks = all_peaks)
saveRDS(other_data, here::here("data", "nsclc_other_peak_data.rds"))
