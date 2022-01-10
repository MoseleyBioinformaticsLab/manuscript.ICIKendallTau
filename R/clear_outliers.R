clear_outliers = function(outlier_data){
  split_method = split(outlier_data, outlier_data[, c("method", "keep_num")])
  
  removed_high = purrr::map_df(split_method, function(in_method){
    max_notout = in_method %>%
      dplyr::filter(!outlier) %>%
      dplyr::pull(med_cor) %>%
      max()
    in_method[in_method$med_cor >= max_notout, "outlier"] = FALSE
    in_method
  })
}
