##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param select_random_fraction
##' @param transcript_pca
##' @return
##' @author rmflight
##' @export
evaluate_by_pca <- function(select_random_fraction, transcript_pca) {

  used_select = rownames(select_random_fraction$data)
  loading_mat = transcript_pca$rotation
  selected_frac = purrr::map_dbl(seq(1, ncol(loading_mat)), function(use_pc){
    pc_total = sum(abs(loading_mat[, use_pc]))
    select_total = sum(abs(loading_mat[used_select, use_pc]))
    frac = select_total / pc_total
    frac
  })
  
  perc_var = visqc_score_contributions(transcript_pca$x)
  perc_var$selected_frac = selected_frac
  perc_var$type = select_random_fraction$type
  perc_var$frac = select_random_fraction$frac
  perc_var
}

summarize_by_pca = function(perc_var){
  perc_var %>%
    dplyr::group_by(type, frac) %>%
    dplyr::summarise(kt = ici_kt(percent, selected_frac)["tau"],
                     med = median(selected_frac),
                     diff = med - frac[1],
                     diff_perc = diff / frac[1])
}
