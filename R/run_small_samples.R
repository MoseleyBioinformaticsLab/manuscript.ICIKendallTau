##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param select_ss_small
##' @return
##' @author rmflight
##' @export
run_small_samples <- function(select_ss_small) {

  future::plan(multicore(workers = 1))
  cor = visqc_ici_kendallt(t(in_data$data))
  cor$type = in_data$type
  cor$n_sample = in_data$n_sample
  cor
  future::plan(multicore)

}
