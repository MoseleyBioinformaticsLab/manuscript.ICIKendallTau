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

right_censor_correlate = function(rc_samples, cut_values = seq(0, 1.5, by = 0.1)){
  censor_cor = purrr::map_df(cut_values, function(in_cut){
    tmp_rc = rc_samples
    tmp_rc[tmp_rc < in_cut] = NA
    n_na = sum(is.na(tmp_rc))
    
    ici_cor = ici_kendalltau(t(tmp_rc), exclude_0 = FALSE, perspective = "global")$cor[1,2]
    p_cor = cor(tmp_rc, use = "pairwise.complete.obs", method = "pearson")[1,2]
    k_cor = cor(tmp_rc, use = "pairwise.complete.obs", method = "kendall")[1,2]
    tmp_rc[is.na(tmp_rc)] = 0
    p_cor_0 = cor(tmp_rc, method = "pearson")[1,2]
    k_cor_0 = cor(tmp_rc, method = "kendall")[1,2]
    
    data.frame(cor = c(ici_cor,
                       p_cor,
                       k_cor,
                       p_cor_0,
                       k_cor_0),
               which = c("ici",
                         "pearson",
                         "kendall",
                         "pearson_0",
                         "kendall_0"),
               n_na = n_na,
               cutoff = in_cut)
    
  })
  censor_cor
}

random_censor_correlate = function(rc_samples, n_na = seq(0, 300, 50)){
  censor_cor = purrr::map_df(n_na, function(in_na){
    tmp_rc = rc_samples
    n_total = length(rc_samples)
    na_locs = sample(n_total, in_na)
    tmp_rc[na_locs] = NA
    
    ici_cor = ici_kendalltau(t(tmp_rc), exclude_0 = FALSE, perspective = "global")$cor[1,2]
    p_cor = cor(tmp_rc, use = "pairwise.complete.obs", method = "pearson")[1,2]
    k_cor = cor(tmp_rc, use = "pairwise.complete.obs", method = "kendall")[1,2]
    tmp_rc[is.na(tmp_rc)] = 0
    p_cor_0 = cor(tmp_rc, method = "pearson")[1,2]
    k_cor_0 = cor(tmp_rc, method = "kendall")[1,2]
    
    tmp_frame = data.frame(cor = c(ici_cor,
                       p_cor,
                       k_cor,
                       p_cor_0,
                       k_cor_0),
               which = c("ici",
                         "pearson",
                         "kendall",
                         "pearson_0",
                         "kendall_0"),
               n_na = in_na,
               cutoff = 0)
    tmp_frame$na_locs = list(na_locs)
    tmp_frame
    
  })
  censor_cor
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
