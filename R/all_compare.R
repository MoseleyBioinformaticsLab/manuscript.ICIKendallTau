##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
all_kendalltau <- function(s1, s2, sneg, where_na, perspective = "global") {
  n_entry = length(s1)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  tmp = furrr::future_map_dfr(where_na, function(use_na){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_s1 = s1
    tmp_s2 = s2
    tmp_neg = sneg
    
    s2_na = use_na[use_na > n_entry] - n_entry
    s1_na = use_na[use_na <= n_entry]
    neg_na = (n_entry + 1) - (use_na[use_na > n_entry] - n_entry)
    
    tmp_s1[s1_na] = NA
    tmp_s2[s2_na] = NA
    tmp_neg[neg_na] = NA
    
    ici_s1 = ici_kt(tmp_s1, tmp_s2, perspective = perspective)[[1]]
    ici_neg = ici_kt(tmp_s1, tmp_neg, perspective = perspective)[[1]]
    data.frame(cor = c(ici_s1, ici_neg),
               comp = c("pos", "neg"))
  })
  tmp
}

319715810
