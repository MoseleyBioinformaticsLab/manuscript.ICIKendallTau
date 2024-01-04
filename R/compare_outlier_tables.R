add_method = function(outlier_df, map_method = c("icikt" = "ICI-Kt",
                                                 "icikt_complete" = "ICI-Kt * Completeness",
                                                 "pearson_base" = "Pearson Base",
                                                 "pearson_base_nozero" = "Pearson No Zeros",
                                                 "pearson_log1p" = "Pearson Log(x + 1)",
                                                 "pearson_log" = "Pearson Log(x)",
                                                 "kt_base" = "Kendall-tau")){
  outlier_df$method = ""
  outlier_df = outlier_df[outlier_df$which %in% names(map_method), ]
  for (imap in names(map_method)) {
    match_loc = outlier_df$which %in% imap
    outlier_df[match_loc, "method"] = map_method[imap]
  }
  outlier_df$method = factor(outlier_df$method, levels = map_method, ordered = TRUE)
  outlier_df
}

perc_to_number = function(in_table){
  in_table %>%
  dplyr::mutate(keep_num = dplyr::case_when(
    keep_num < 1 ~ paste0(keep_num * 100, "%"),
    TRUE ~ as.character(keep_num)
  ))
}

compare_outlier_tables = function(outlier_list, keep_compare, sort_var, map_method = c("icikt" = "ICI-Kt",
                                                                                       "icikt_complete" = "ICI-Kt * Completeness",
                                                                                       "pearson_base" = "Pearson Base",
                                                                                       "pearson_base_nozero" = "Pearson No Zeros",
                                                                                       "pearson_log1p" = "Pearson Log(x + 1)",
                                                                                       "pearson_log" = "Pearson Log(x)",
                                                                                       "kt_base" = "Kendall-tau")){
  
  # 
  # outlier_list = yeast_outliers_1$outliers
  # sort_var = "pearson_log"
  # keep_compare = c("icikt", "icikt_complete", "pearson_log")
  # outlier_list = yeast_single2
  # keep_compare = c("icikt", "icikt_complete", "pearson_log", "pearson_base_nozero", "manuscript")
  # sort_var = "manuscript"
  # map_method = c("icikt" = "ICI-Kt",
  #                                                                                                                                                         "icikt_complete" = "ICI-Kt * Completeness",
  #                                                                                                                                                         "pearson_base" = "Pearson Base",
  #                                                                                                                                                         "pearson_base_nozero" = "Pearson No Zeros",
  #                                                                                                                                                         "pearson_log1p" = "Pearson Log(x + 1)",
  #                                                                                                                                                         "pearson_log" = "Pearson Log(x)",
  #                                                                                                                                                         "kt_base" = "Kendall-tau",
  #                                                                                                                                                         "manuscript" = "Manuscript")
  
  outlier_list = add_method(outlier_list, map_method = map_method)
  outlier_list = outlier_list %>%
    dplyr::filter(which %in% keep_compare)
  
  isoutlier = outlier_list %>%
    dplyr::filter(outlier)
  
  if (nrow(isoutlier) == 0) {
    return(NULL)
  }
  
  out_samples = data.frame(sample_id = unique(isoutlier$sample_id))
  
  # reordering by the one we want
  order_tmp = outlier_list %>%
    dplyr::filter(which %in% sort_var) %>%
    dplyr::filter(sample_id %in% out_samples$sample_id)
  order_split = split(order_tmp, order_tmp$sample_class)
  order_isorder = purrr::map_df(order_split, function(in_split){
    in_split %>%
      dplyr::arrange(dplyr::desc(med_cor))
  })
  
  # creating the main table
  out_by_which = split(outlier_list, outlier_list$which)
  out_by_which = out_by_which[keep_compare]
  out_table = purrr::map_dfc(out_by_which, function(in_which){
    tmp_table = dplyr::left_join(out_samples, in_which, by = "sample_id")
    tmp_table = tmp_table %>%
      dplyr::mutate(outlier2 = dplyr::case_when(
        outlier ~ "X",
        !outlier ~ " "
      ))
    new_frame = data.frame(outlier = tmp_table$outlier2,
                           cor = tmp_table$med_cor)
    names(new_frame) = paste0(in_which$which[1], ".", names(new_frame))
    
    new_frame
  })
  out_table$sample_id = out_samples$sample_id
  out_table = dplyr::left_join(order_isorder[, c("sample_id"), drop = FALSE], out_table, by = "sample_id")
  
  ft_out = flextable(out_table)
  new_labels = purrr::map(names(out_table), function(in_name){
    if (grepl("cor", in_name)) {
      return("Correlation")
    }
    if (grepl("sample", in_name)) {
      return("Sample")
    }
    if (grepl("outlier", in_name)) {
      return("Outlier")
    }
  })
  names(new_labels) = names(out_table)
  
  which_method = gsub("[.].*", "", names(out_table))
  which_map = map_method[which_method]
  names(which_map) = NULL
  which_map[1] = ""
  
  rle_which = rle(which_map)
  
  ft_out = set_header_labels(ft_out,
                              values = new_labels)
  ft_out = add_header_row(ft_out,
                            values = rle_which$values,
                            colwidths = rle_which$lengths)
  ft_out = colformat_double(ft_out, digits = 3)
  ft_out
}
