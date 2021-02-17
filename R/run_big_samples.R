##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param in_data
##' @return
##' @author rmflight
##' @export
run_big_samples <- function(in_data) {

  cor = visqc_ici_kendallt(t(in_data$data))
  cor$type = in_data$type
  cor$n_sample = in_data$n_sample
  cor

}
