calculate_qratio = function(network, annotations)
{
  # network = network_correlations
  # annotations = use_annotation
  # assume the network is a data.frame edges
  # sum over annotations
  # sum the entire network of weights
  
  # assume only the positive partial correlations are useful
  network = network |>
    dplyr::filter(weight > 0)
  all_features_annotations = unique(unlist(annotations))
  
  network_all_nodes = unique(unlist(network[, c("start_node", "end_node")]))
  
  network_annotated = network |>
    dplyr::filter(start_node %in% all_features_annotations,
                  end_node %in% all_features_annotations)
  
  # network_annotated = network_annotated |>
  #   dplyr::mutate(weight = 1)
  
  annotated_sum = sum(network_annotated$weight)
  partitions = purrr::imap(annotations, \(in_features, id){
    other_features = all_features_annotations[!(all_features_annotations %in% in_features)]
    # generate self-self network - sum correlations
    within_sum = network_annotated |>
      dplyr::filter(((start_node %in% in_features) & (end_node %in% in_features)) 
                    | ((end_node %in% in_features) & (start_node %in% in_features))) |>
      dplyr::summarise(sum = sum(weight)) |>
      dplyr::pull(sum)
    out_sum = network_annotated |>
      dplyr::filter(((start_node %in% in_features) & (end_node %in% other_features)) |
                      (end_node %in% in_features) & (start_node %in% other_features)) |>
      dplyr::summarise(sum = sum(weight)) |>
      dplyr::pull(sum)
    
    out_val = (within_sum / annotated_sum) - ((out_sum / annotated_sum)^2)
    
    n_annot = length(in_features)
    
    tibble::tibble(id = id, ratio = out_val, n_features = length(in_features))
  }) |>
    dplyr::bind_rows()
  q_value = sum(partitions$ratio)
  list(partitions = partitions, q_value = q_value)
}

calculate_feature_network_qratio = function(partial_correlations, annotations, 
                                            compound_type = "pathway",
                                            compound_mapping = metabolite_kegg)
{
  # partial_correlations = tar_read(feature_partial_cor_ici_yeast)
  # partial_correlations = tar_read(feature_partial_cor_ici_ratstamina)
  # annotations = tar_read(feature_annotations)
  # compound_type = "pathway"
  # compound_mapping = tar_read(metabolite_kegg)
  # 
  if (is.na(partial_correlations$pcor[["partial_cor"]][1])) {
    network_partitioning = list(partitions = tibble::tibble(id = NA, ratio = NA, n_features = NA),
                                q_value = NA)
  } else {
    
    network_correlations = partial_correlations$pcor |>
      dplyr::filter(significant) |>
      dplyr::transmute(start_node = s1,
                       end_node = s2,
                       weight = partial_cor)
    
    if (grepl("yeast", partial_correlations$data_id)) {
      use_annotation = annotations$yeast@annotation_features
    } else if (grepl("rat", partial_correlations$data_id)) {
      
      network_correlations2 = dplyr::left_join(
        network_correlations, compound_mapping[, c("feature_id", "kegg_id2")], 
        by = c("start_node" = "feature_id"),
        relationship = "many-to-many"
      ) |>
        dplyr::mutate(s_kegg = kegg_id2,
                      kegg_id2 = NULL)
      network_correlations2 = dplyr::left_join(
        network_correlations2, compound_mapping[, c("feature_id", "kegg_id2")],
        by = c("end_node" = "feature_id"),
        relationship = "many-to-many"
      ) |> 
        dplyr::mutate(e_kegg = kegg_id2)
      network_correlations = network_correlations2 |>
        dplyr::transmute(start_node = s_kegg,
                         end_node = e_kegg,
                         weight = weight)
      
      if (compound_type %in% "pathway") {
        use_annotation = annotations$kegg_pathway@annotation_features
        
      }
    } else if (grepl("adenocarcinoma", partial_correlations$data_id)) {
      use_annotation = annotations$human@annotation_features
    } else if (grepl("egfrgenotype", partial_correlations$data_id)) {
      use_annotation = annotations$mouse@annotation_features
    } else if (grepl("nsclc", partial_correlations$data_id)) {
      use_annotation = annotations$lipid_class@annotation_features
    }
    
    network_partitioning = calculate_qratio(network_correlations, use_annotation)
  }
  
  
  network_partitioning$data_id = partial_correlations$data_id
  network_partitioning$method_id = partial_correlations$method_id
  network_partitioning$full_id = partial_correlations$full_id
  return(network_partitioning)
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
  # feature_data = tar_read(feature_correlation_pearson_base_nozero_ratstamina)
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
  pcor_vals = try(cor_to_pcor(cor_matrix))
  
  if (inherits(pcor_vals, "try-error")) {
    pcor_p_values = tibble::tibble(partial_cor = NA, id = NA)
  } else {
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
  }
  
  #pcor_vals[upper.tri(pcor_vals)] = NA
  
  return(list(pcor = pcor_p_values,
       data_id = feature_data$data_id,
       method_id = feature_data$method_id,
       full_id = feature_data$full_id))
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

get_feature_annotations = function(kegg_data)
{
  compound_pathway = compound_annotation(kegg_data, "pathway")
  compound_network = compound_annotation(kegg_data, "network")
  compound_module = compound_annotation(kegg_data, "module")
  
  yeast_keys = AnnotationDbi::keys(org.Sc.sgd.db)
  yeast_entrez = suppressMessages(AnnotationDbi::select(org.Sc.sgd.db, keys = yeast_keys, columns = "ENTREZID"))
  yeast_uniq_entrez = unique(yeast_entrez$ENTREZID)
  yeast_uniq_entrez = yeast_uniq_entrez[!is.na(yeast_uniq_entrez)]
  yeast_reactome = suppressMessages(AnnotationDbi::select(reactome.db, keys = yeast_uniq_entrez,
                                                          keytype = "ENTREZID", columns = c("PATHNAME", "REACTOMEID")))
  yeast_entrez = yeast_entrez |>
    dplyr::transmute(ID = ORF,
                     ENTREZID = ENTREZID)
  yeast_id_reactome = dplyr::left_join(yeast_entrez, yeast_reactome, by = "ENTREZID")
  
  yeast_annotation = create_feature_annotation_object(yeast_id_reactome)
  
  mouse_keys = AnnotationDbi::keys(org.Mm.eg.db)
  mouse_entrez = suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = mouse_keys, columns = "ENSEMBL"))
  mouse_uniq_entrez = unique(mouse_entrez$ENTREZID)
  mouse_uniq_entrez = mouse_uniq_entrez[!is.na(mouse_uniq_entrez)]
  mouse_reactome = suppressMessages(AnnotationDbi::select(reactome.db, keys = mouse_uniq_entrez,
                                                          keytype = "ENTREZID", columns = c("PATHNAME", "REACTOMEID")))
  mouse_entrez = mouse_entrez |>
    dplyr::transmute(ID = ENSEMBL,
                     ENTREZID = ENTREZID)
  mouse_id_reactome = dplyr::left_join(mouse_entrez, mouse_reactome, by = "ENTREZID", relationship = "many-to-many")
  mouse_annotation = create_feature_annotation_object(mouse_id_reactome)
  
  human_keys = AnnotationDbi::keys(org.Hs.eg.db)
  human_entrez = suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = human_keys, columns = "ENSEMBL"))
  human_uniq_entrez = unique(human_entrez$ENTREZID)
  human_uniq_entrez = human_uniq_entrez[!is.na(human_uniq_entrez)]
  human_reactome = suppressMessages(AnnotationDbi::select(reactome.db, keys = human_uniq_entrez,
                                                          keytype = "ENTREZID", columns = c("PATHNAME", "REACTOMEID")))
  human_entrez = human_entrez |>
    dplyr::transmute(ID = ENSEMBL,
                     ENTREZID = ENTREZID)
  human_id_reactome = dplyr::left_join(human_entrez, human_reactome, by = "ENTREZID", relationship = "many-to-many")
  human_annotation = create_feature_annotation_object(human_id_reactome)
  
  list(human = human_annotation,
       mouse = mouse_annotation,
       yeast = yeast_annotation,
       kegg_pathway = compound_pathway,
       kegg_module = compound_module,
       kegg_network = compound_network)
}

create_feature_annotation_object = function(annotation_df)
{
  # annotation_df = yeast_id_reactome
  annotation_df = annotation_df |>
    dplyr::filter(!is.na(REACTOMEID))
  split_data = split(annotation_df$ID, annotation_df$REACTOMEID) |>
    purrr::map(unique)
  
  pathway_meta = annotation_df |>
    dplyr::select(PATHNAME, REACTOMEID) |>
    dplyr::distinct()
  pathway_description = pathway_meta[["PATHNAME"]]
  names(pathway_description) = pathway_meta[["REACTOMEID"]]
  pathway_description = pathway_description[names(split_data)]
  annotation_obj = categoryCompare2::annotation(split_data,
                                                description = pathway_description)
  annotation_obj
}
