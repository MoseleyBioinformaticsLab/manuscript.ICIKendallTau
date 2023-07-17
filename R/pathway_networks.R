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

create_networks = function(feature_correlations)
{
  # feature_correlations = tar_read(yeast_feature_cor)
  # assume correlations have been calculated all the different ways
  various_network = furrr::future_map(feature_correlations, \(in_cor){
    # in_cor = feature_correlations$pearson_log1p
    # in_cor = in_cor[1:1000, 1:1000]
    # some of our kendall-tau correlations do not have 1 on the diagonal
    diag(in_cor) = 1
    pcor_vals = cor_to_pcor(in_cor)
    network = tidygraph:::as_tbl_graph.matrix(pcor_vals, directed = FALSE, diag = FALSE)
    all_edges = network |>
      activate(edges) |>
      as_tibble() |>
      dplyr::mutate(padj = p.adjust(1 - abs(weight), method = "bonferroni"))
    
  })
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
