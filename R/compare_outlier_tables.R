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
  outlier_df
}

compare_outlier_tables = function(outlier_list, keep_compare, sort_var, map_method = c("icikt" = "ICI-Kt",
                                                                                       "icikt_complete" = "ICI-Kt * Completeness",
                                                                                       "pearson_base" = "Pearson Base",
                                                                                       "pearson_base_nozero" = "Pearson No Zeros",
                                                                                       "pearson_log1p" = "Pearson Log(x + 1)",
                                                                                       "pearson_log" = "Pearson Log(x)",
                                                                                       "kt_base" = "Kendall-tau")){
  
  # drake::loadd(yeast_outliers_1)
  # outlier_list = yeast_outliers_1$outliers
  # sort_var = "pearson_log"
  # keep_compare = c("icikt", "icikt_complete", "pearson_log")
  # 
  
  outlier_list = add_method(outlier_list, map_method = map_method)
  outlier_list = outlier_list %>%
    dplyr::filter(which %in% keep_compare)
  
  isoutlier = outlier_list %>%
    dplyr::filter(outlier)
  
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
  out_by_method = split(outlier_list, outlier_list$method)
  out_table = purrr::map_dfc(out_by_method, function(in_method){
    tmp_table = dplyr::left_join(out_samples, in_method, by = "sample_id")
    tmp_table = tmp_table %>%
      dplyr::mutate(outlier2 = dplyr::case_when(
        outlier ~ "X",
        !outlier ~ " "
      ))
    new_frame = data.frame(outlier = tmp_table$outlier2,
                           cor = tmp_table$med_cor)
    names(new_frame) = paste0(in_method$which[1], ".", names(new_frame))
    
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
  
  ft_out = set_header_labels(ft_out,
                              values = new_labels)
  header_list = map_method[keep_compare]
  names(header_list) = NULL
  header_list = c("", header_list)
  ft_out = add_header_row(ft_out,
                            values = header_list,
                            colwidths = c(1, rep(2, length(mapping_list))))
  ft_out = colformat_double(ft_out, digits = 2)
  list(plot_data = outliers_use,
       table = ft_out)
}
