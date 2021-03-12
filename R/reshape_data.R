##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param combined_random
##' @return
##' @author rmflight
##' @export
reshape_data <- function(combined_data) {
  # this line is for testing only!!
  #combined_data = readRDS("data/random_results.rds")
  combined_data = combined_data %>%
    dplyr::mutate(logtime = log(time),
                  logsd = -1 * log10(sd))
  combined_long = combined_data %>%
    tidyr::pivot_longer(!c(fraction, type), names_to = "measure", values_to = "value")
  combined_long
}
