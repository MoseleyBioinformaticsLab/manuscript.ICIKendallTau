##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author rmflight
##' @export
create_rc_samples <- function() {

  base_sample = rlnorm(1000, meanlog = 1, sdlog = 0.5)
  rep_data = add_uniform_noise(2, base_sample, 0.2)
  rep_data

}

add_uniform_noise <- function(n_rep, value, sd, use_zero = FALSE){
  n_value <- length(value)
  
  n_sd <- n_rep * n_value
  
  out_sd <- rnorm(n_sd, 0, sd)
  out_sd <- matrix(out_sd, nrow = n_value, ncol = n_rep)
  
  if (!use_zero){
    tmp_value <- matrix(value, nrow = n_value, ncol = n_rep, byrow = FALSE)
    out_value <- tmp_value + out_sd
  } else {
    out_value <- out_sd
  }
  
  return(out_value)
}
