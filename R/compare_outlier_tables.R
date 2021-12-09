compare_outlier_tables = function(outlier_list, mapping_list, sort_var){
  # drake::loadd(yeast_outliers)
  # outlier_list = yeast_outliers
  # sort_var = "pearson_log"
  # mapping_list = c("icikt" = "ICI-Kt", 
                   # "icikt_complete" = "ICI-Kt * Completeness",
                   # "pearson_log" = "Pearson")
  outliers_use = purrr::map_df(names(mapping_list), function(in_id){
    #message(in_id)
    tmp_out = outlier_list[[in_id]]
    tmp_out$method = mapping_list[[in_id]]
    tmp_out
  })
  
  isoutlier = outliers_use %>%
    dplyr::filter(outlier)
  
  out_samples = data.frame(sample_id = unique(isoutlier$sample_id))
  
  # reordering by the one we want
  order_tmp = outlier_list[[sort_var]] %>%
    dplyr::filter(sample_id %in% out_samples$sample_id)
  order_split = split(order_tmp, order_tmp$sample_class)
  order_isorder = purrr::map_df(order_split, function(in_split){
    in_split %>%
      dplyr::arrange(dplyr::desc(med_cor))
  })
  
  # creating the main table
  out_by_method = split(outliers_use, outliers_use$method)
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
  header_list = mapping_list
  names(header_list) = NULL
  header_list = c("", header_list)
  ft_out = add_header_row(ft_out,
                            values = header_list,
                            colwidths = c(1, rep(2, length(mapping_list))))
  ft_out = colformat_double(ft_out, digits = 2)
  list(plot_data = outliers_use,
       table = ft_out)
}
