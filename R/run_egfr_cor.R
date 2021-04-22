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
