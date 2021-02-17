##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param transcript_data
##' @param n_sample
##' @return
##' @author rmflight
##' @export
select_samples <- function(transcript_data, n_sample = n_sample) {
  use_cols = sample(ncol(transcript_data), n_sample)
  
  list(data = transcript_data[, use_cols],
       type = "random",
       n_sample = n_sample)

}
