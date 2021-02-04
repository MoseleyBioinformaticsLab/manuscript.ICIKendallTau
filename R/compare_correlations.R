##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_positive_kt <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  tmp = furrr::future_map_dbl(where_na, function(use_na){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    ici_kendallt(tmp_x, tmp_y, perspective = perspective)
  }, .progress = TRUE)
  
  tmp
}

compare_positive_kt_c <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  tmp = furrr::future_map_dbl(where_na, function(use_na){
  #purrr::imap_dbl(where_na, function(use_na, .y){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    ici_kendallt(tmp_x, tmp_y, perspective = perspective)
  #})
  }, .progress = TRUE)
  
  tmp
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_negative_kt <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  furrr::future_map_dbl(where_na, function(use_na){
    #message(.y)
    tmp_x = x
    tmp_y = y
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    y_na = n_entry - y_na + 1
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    ici_kendallt(tmp_x, tmp_y, perspective = perspective)
  }, .progress = TRUE)
  
  
}

compare_negative_kt_c <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  furrr::future_map_dbl(where_na, function(use_na){
    #message(.y)
    tmp_x = x
    tmp_y = y
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    y_na = n_entry - y_na + 1
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    ici_kendallt(tmp_x, tmp_y, perspective = perspective)
  }, .progress = TRUE)
  
  
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_negative_pearson <- function(x, y, where_na, low_indices = FALSE) {
  n_entry = length(x)
  furrr::future_map_dbl(where_na, function(use_na){
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    y_na = n_entry - y_na + 1
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    in_matrix = rbind(tmp_x, tmp_y)
    out_res = locally_it_weighted_pairwise_correlation(in_matrix)
    out_res$cor[1, 2]
  }, .progress = TRUE)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_positive_pearson <- function(x, y, where_na, low_indices = FALSE) {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  furrr::future_map_dbl(where_na, function(use_na){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    in_matrix = rbind(tmp_x, tmp_y)
    out_res = locally_it_weighted_pairwise_correlation(in_matrix)
    out_res$cor[1, 2]
  }, .progress = TRUE)
  
  
}
