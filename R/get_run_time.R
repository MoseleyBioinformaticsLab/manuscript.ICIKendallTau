##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param run_small
##' @return
##' @author rmflight
##' @export
get_run_time <- function(in_cor) {

  data.frame(time = in_cor$run_time,
             n_sample = in_cor$n_sample)

}
