process_study_data = function(in_data)
{
  # in_data = tar_read(study_data)[[4]]
  subject_data = in_data$SUBJECT_SAMPLE_FACTORS |>
    janitor::clean_names()
  subject_more = purrr::imap(subject_data, function(in_data, in_name){
    if (inherits(in_data, "data.frame")) {
      tmp_frame = in_data |>
        janitor::clean_names()
      names(tmp_frame) = paste0(in_name, ".", names(tmp_frame))
      return(tmp_frame)
    } else {
      return(in_data)
    }
  }) |>
    dplyr::bind_cols()
  subject_more$sample = janitor::make_clean_names(subject_more$sample_id)
  
  factor_cols = which(grepl("^factor", colnames(subject_more)))
  
  subject_more$treatment = paste0(subject_more[[factor_cols[1]]], ":", subject_more[[factor_cols[2]]])
  
  measured_name = which(grepl("METABOLITE_DATA", names(in_data)))
  if (length(measured_name) == 0) {
    return(NULL)
  }
  measurement_data = in_data[[measured_name]][["Data"]] |>
    janitor::clean_names()
  measurement_only = measurement_data |>
    dplyr::select(-metabolite)
  measurement_values = purrr::map(measurement_only, function(in_data){
    tmp = as.numeric(in_data)
    tmp[is.na(tmp)] = 0
    tmp
  }) |>
    dplyr::bind_cols()
  normed_values = normalize_median(measurement_values) |>
    as.matrix()
  tmp_features = gsub("\\^|`", "'", measurement_data$metabolite)
  rownames(normed_values) = tmp_features
  
  # prefer kegg ID if it exists, and map only the **first** one if there are multiple
  alias_data = in_data$MS_METABOLITE_DATA$Metabolites |>
    janitor::clean_names() |>
    dplyr::select(metabolite, kegg_id) |>
    dplyr::transmute(feature_id = gsub("\\^|`", "'", metabolite),
                     trimmed_id = gsub("_", ",", gsub("_1$|_2$", "", feature_id)),
                     kegg_id = kegg_id)
  
  alias_split = split(alias_data, alias_data$trimmed_id)
  
  alias_unique = purrr::map(alias_split, \(in_split){
    if (nrow(in_split) == 1) {
      return(in_split)
    }
    any_kegg_zero = nchar(in_split$kegg_id) == 0
    if (any(any_kegg_zero)) {
      non_zero = in_split$kegg_id[!any_kegg_zero]
      if (length(unique(non_zero)) == 1) {
        use_kegg = unique(non_zero)
        in_split[["kegg_id"]][any_kegg_zero] = use_kegg
      } else {
        in_split[["kegg_id"]][any_kegg_zero] = in_split[["trimmed_id"]][any_kegg_zero]
      }
      return(in_split |> dplyr::distinct())
    } else {
      return(in_split |> dplyr::distinct())
    }
  }) |>
    dplyr::bind_rows()
  
  alias_unique2 = alias_unique |>
    dplyr::select(trimmed_id, kegg_id) |>
    dplyr::distinct()
  
  row_aliases = tibble::tibble(feature_id = rownames(normed_values),
                               trimmed_id = gsub("_", ",", gsub("_1$|_2$", "", feature_id)))
  
  row_aliases_kegg = dplyr::left_join(row_aliases, alias_unique2, by = "trimmed_id")
  row_aliases_kegg = row_aliases_kegg |>
    dplyr::mutate(kegg_id2 = dplyr::case_when(
      nchar(kegg_id) == 0 ~ trimmed_id,
      nchar(kegg_id) > 0 ~ kegg_id
    ))
  
  return(list(counts = normed_values,
              info = subject_more,
              metabolite_mapping = row_aliases_kegg))
}


normalize_median = function(sample_measure_values)
{
  out_norm = purrr::imap(sample_measure_values, function(in_values, col_id){
    tmp_values = in_values
    tmp_values = tmp_values[in_values > 0]
    median_value = median(tmp_values)
    out_values = in_values / median_value
    out_values
  }) |>
    dplyr::bind_cols()
  return(out_norm)
}

pos_data = jsonlite::fromJSON(here::here("raw_data", "mwtab", "ST000017_AN000034.json"))
pos_processed_data = process_study_data(pos_data)
pos_processed_data$data_id = "mwtab_ratstamina"
pos_metabolite_mapping = pos_processed_data$metabolite_mapping
pos_processed_data$metabolite_mapping = NULL

saveRDS(pos_processed_data, here::here("data/mwtab_st000017_an000034_count_info.rds"))
saveRDS(pos_metabolite_mapping, here::here("data/mwtab_st000017_an000034_metabolites.rds"))
