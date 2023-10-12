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

pcor_pvalue_z = function(pcor_values)
{
  mean_pcor = mean(pcor_values)
  sd_pcor = sd(pcor_values)
  z_scores = (pcor_values - mean_pcor) / sd_pcor
  
  p_value = 2 * pnorm(-abs(z_scores))
  p_adjust = p.adjust(p_value, method = "BH")
  tibble::tibble(pvalue = p_value,
                 padjust = p_adjust)
}

pcor_pvalue_extreme = function(pcor_value_df, p_cut = 0.05)
{
  # pcor_value_df = pcor_long
  # p_cut = 0.05
  mean_partial = mean(pcor_value_df$partial_cor)
  pcor_value_df = pcor_value_df |>
    dplyr::arrange(partial_cor) |>
    dplyr::mutate(rank = rank(partial_cor),
                  fraction = rank / max(rank),
                  p_value = dplyr::case_when(
                    partial_cor < mean_partial ~ fraction,
                    partial_cor > mean_partial ~ 1 - fraction
                  ),
                  significant = p_value <= (p_cut / 2))
  pcor_value_df
}

n_extreme = function(in_value_df)
{
  in_value_df = in_value_df |>
    dplyr::mutate(rank = rank(partial_cor),
                  p_value = 1 - ((max(rank) - rank) / max(rank)))
}


calculate_partial_cor_pvalues = function(feature_data)
{
  # feature_data = tar_read(feature_correlation_ici_yeast)
  # feature_data = tar_read(feature_correlation_pearson_base_ratstamina)
  # feature_data = tar_read(feature_correlation_ici_completeness_ratstamina)
  if (inherits(feature_data$cor, "data.frame")) {
    feature_correlations = feature_data$cor |>
      dplyr::transmute(s1 = s1,
                       s2 = s2,
                       cor = raw,
                       id = paste0(s1, ".", s2))
  } else {
    feature_correlations = feature_data$cor$cor
    feature_correlations = feature_correlations |>
      dplyr::mutate(id = paste0(s1, ".", s2)) |>
      dplyr::arrange(id)
    if (!is.null(feature_data$completeness)) {
      completeness = feature_data$completeness |>
        dplyr::mutate(id = paste0(s1, ".", s2)) |>
        dplyr::arrange(id)
      if (all.equal(completeness$id, feature_correlations$id)) {
        feature_correlations$cor = feature_correlations$cor * completeness$completeness
      }
    }
  }
  
  
  message("converting to matrix form ...\n")
  cor_matrix = long_df_2_cor_matrix(feature_correlations |> dplyr::select(s1, s2, cor))
  diag(cor_matrix) = 1
  cor_matrix[is.na(cor_matrix)] = 0
  message("calculating partial correlation ...\n")
  pcor_vals = cor_to_pcor(cor_matrix)
  
  #pcor_vals[upper.tri(pcor_vals)] = NA
  pcor_long = cor_matrix_2_long_df(pcor_vals)
  pcor_long = pcor_long |>
    dplyr::filter(!is.na(cor)) |>
    dplyr::mutate(partial_cor = cor,
                  cor = NULL,
                  id = paste0(s1, ".", s2))
  pcor_long = dplyr::left_join(feature_correlations[, c("cor", "id")], pcor_long, by = "id")
  pcor_long = pcor_long |>
    dplyr::filter(!(s1 == s2))
  message("adjusting p-values")
  pcor_long = pcor_long |>
    dplyr::filter(!is.na(partial_cor))
  pcor_p_values = pcor_pvalue_extreme(pcor_long)
  list(pcor = pcor_p_values,
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

compound_annotation = function(data_file, type = "pathway")
{
  # data_file = "data/kegg_compound_mapping.rds"
  # type = "pathway"
  compound_data = readRDS(data_file)
  
  use_entries = compound_data[[type]]
  split_data = split(use_entries$compound, use_entries$id) |>
    purrr::map(unique)
  
  pathway_meta = compound_data$meta |>
    dplyr::filter(type %in% type)
  pathway_description = pathway_meta[["description"]]
  names(pathway_description) = pathway_meta$id
  pathway_description = pathway_description[names(split_data)]
  annotation_obj = categoryCompare2::annotation(split_data, annotation_type = paste0("kegg-", type),
                                                description = pathway_description,
                                                feature_type = "kegg-compound")
  annotation_obj
}
