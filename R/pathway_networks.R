calculate_qratio = function(network, annotations)
{
  # assume the network is a data.frame of nodes and edges
  # sum over annotations
  # sum the entire network of weights
  # network_sum = sum(network$edges$weights)
  p_vals = purrr::map_dbl(annotations, \(in_features){
    # generate self-self network - sum correlations
    in_network = network |>
      dplyr::filter(start_node %in% in_features, end_node %in% in_features)
    in_sum = sum(in_network$edges$weights)
    out_network = network |>
      dplyr::filter(start_node %in% in_features)
    out_sum = sum(in_network$edges$weights)
    
    out_val = (in_sum / network_sum) - ((out_sum / network_sum)^2)
    
    out_val
  })
  q = sum(p_vals)
  q
}

pcor_pvalue = function(pcor_values)
{
  mean_pcor = mean(pcor_values)
  sd_pcor = sd(pcor_values)
  z_scores = (pcor_values - mean_pcor) / sd_pcor
  
  p_value = 2 * pnorm(-abs(z_scores))
  p_adjust = p.adjust(p_value, method = "BH")
  tibble::tibble(pvalue = p_value,
                 padjust = p_adjust)
}


calculate_pcor_pvalues = function(feature_data)
{
  # feature_data = tar_read(feature_correlation_ici_yeast)
  if (inherits(feature_data$cor, "data.frame")) {
    feature_correlations = feature_data$cor
  } else {
    feature_correlations = feature_data$cor$cor
  }
  
  message("converting to matrix form ...\n")
  cor_matrix = long_df_2_cor_matrix(feature_correlations |> dplyr::select(s1, s2, cor))
  diag(cor_matrix) = 1
  cor_matrix[is.na(cor_matrix)] = 0
  message("calculating partial correlation ...\n")
  pcor_vals = cor_to_pcor(cor_matrix)
  
  pcor_vals[upper.tri(pcor_vals, diag = TRUE)] = NA
  pcor_long = cor_matrix_2_long_df(pcor_vals)
  pcor_long = pcor_long |>
    dplyr::mutate(partial_cor = cor,
                  cor = NULL)
  message("adjusting p-values")
  pcor_long = pcor_long |>
    dplyr::filter(!is.na(partial_cor))
  pcor_long = dplyr::bind_cols(pcor_long,
                               pcor_pvalue(pcor_long$partial_cor))
  list(pcor = pcor_long,
       data_id = feature_data$data_id,
       method_id = feature_data$method_id,
       full_id = feature_data$full_id)
}

create_network_from_pcor_pvalues = function(pcor_data, padjust_limit = 0.05)
{
  # pcor_data = tar_read(feature_pcor_ici_yeast)
  # padjust_limit = 0.05
  pcor_values = pcor_data$pcor
  
  pcor_keep = pcor_values |>
    dplyr::filter(padjust <= padjust_limit)
  
  network = tidygraph::as_tbl_graph(pcor_keep, directed = FALSE)
    all_weights = network |>
      activate(edges) |>
      as_tibble() |>
      dplyr::pull(weight)
    
    p_values = pcor_pvalue(all_weights)
    
    network = network |>
      activate(edges) |>
      cbind(p_values)
    
    mean_weight = mean(all_weights)
    
  
  various_network
}



create_correlation_networks = function(counts_info, keep_num, sample_col, class_col){
  # counts_info = tar_read(yeast_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  counts_completeness = pairwise_completeness(counts_filter)
  counts_cor = run_cor_everyway(t(counts_filter), counts_completeness)
  counts_cor
}

cpp_wrapper = function(code) {
  out <- Rcpp::sourceCpp(code)
  body(out) <- rlang::call2("{", code, body(out))
  out
}

