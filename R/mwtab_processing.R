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
  subject_more$sample_id = janitor::make_clean_names(subject_more$sample_id)
  
  factor_cols = which(grepl("^factor", colnames(subject_more)))
  
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
  normed_values = normalize_median(measurement_values)
  measurement_values$metabolite = measurement_data$metabolite
  
  normed_values$metabolite = measurement_data$metabolite
  
  if (any(grepl("MS", names(in_data)))) {
    analytical_method = "MS"
  } else if (any(grepl("NMR", names(in_data)))) {
    analytical_method = "NMR"
  } else {
    analytical_method = "UNKNOWN"
  }
  return(list(analytical_method = analytical_method,
              project = in_data$PROJECT,
              sample_info = subject_more,
              raw_data = measurement_values,
              normalized_data = normed_values))
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

