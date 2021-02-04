var_select = function(pca_data, matrix_data, fraction){
  n_item = round(fraction * nrow(matrix_data))
  var_rows = apply(matrix_data, 1, var)
  var_data = data.frame(row_id = rownames(matrix_data),
                        var = var_rows)
  var_data = dplyr::arrange(var_data, dplyr::desc(var))
  use_rows = var_data$row_id[1:n_item]
  list(data = matrix_data[use_rows, ],
       type = "var",
       frac = fraction)
}

pca_select = function(pca_data, matrix_data, fraction){
  n_item = round(fraction * nrow(matrix_data))
  pc_cont = visqc_score_contributions(as.matrix(pca_data$x))
  use_pcs = dplyr::filter(pc_cont, cumulative <= 0.95)
  use_pcs$n_row = round(n_item * use_pcs$percent)
  use_rows = vector("character", n_item)
  
  if (is.null(rownames(matrix_data))) {
    tmp_names = seq(1, nrow(matrix_data))
  } else {
    tmp_names = rownames(matrix_data)
  }
  
  count_load = 1
  use_pc = 1
  while (count_load < n_item) {
    use_loadings = pca_data$rotation[, use_pc]
    names(use_loadings) = tmp_names
    use_loadings = sort(abs(use_loadings), decreasing = TRUE)
    use_loadings = use_loadings[!(names(use_loadings) %in% use_rows)]
    use_loc = seq(count_load, count_load + use_pcs$n_row[use_pc] - 1)
    use_loc = use_loc[use_loc <= n_item]
    load_loc = seq(1, length(use_loc))
    use_rows[use_loc] = names(use_loadings[load_loc])
    count_load = use_loc[length(use_loc)] + 1
  }
  list(data = matrix_data[use_rows, ],
       type = "pca",
       frac = fraction)
}

run_fractional_correlation = function(in_data){
  cor = visqc_ici_kendallt(t(in_data$data))
  cor$type = in_data$type
  cor$frac = in_data$frac
  cor
}
