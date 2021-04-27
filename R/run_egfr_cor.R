##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author rmflight
##' @export
run_egfr_cor <- function() {

  count_matrix = readRDS(here::here("data", "brainson_egfr_counts.rds"))
  count_cor = visqc_ici_kendallt(t(count_matrix), exclude_0 = TRUE, perspective = "global")$cor
  count_cor

}

find_egfr_outliers = function(egfr_cor){
  count_matrix = readRDS(here::here("data", "brainson_egfr_counts.rds"))
  egfr_info = readRDS(here::here("data", "brainson_egfr_info.rds"))
  count_matrix = count_matrix[, egfr_info$sample]
  
  egfr_cor = egfr_cor[egfr_info$sample, egfr_info$sample]
  med_cor = median_correlations(egfr_cor, egfr_info$type)
  out_frac = outlier_fraction(t(log1p(count_matrix)), egfr_info$type)
  out_samples = determine_outliers(med_cor, out_frac, frac_weight = 0)
  hi_cor = dplyr::filter(out_samples, med_cor >= median(med_cor), outlier) %>%
    dplyr::pull(sample_id)
  out_samples[out_samples$sample_id %in% hi_cor, "outlier"] = FALSE
  out_samples
}
