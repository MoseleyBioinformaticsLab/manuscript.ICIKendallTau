var_select = function(pca_data, matrix_data, fraction){
  n_item = round(fraction * nrow(matrix_data))
  var_rows = apply(matrix_data, 1, var, na.rm = TRUE)
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

select_random = function(matrix_data, fraction){
  n_item = round(fraction * nrow(matrix_data))
  use_rows = sample(nrow(matrix_data), n_item)
  list(data = matrix_data[use_rows, ],
       type = "random",
       frac = fraction)
}

run_fractional_correlation = function(in_data){
  cor = ici_kendalltau(t(in_data$data))
  cor$type = in_data$type
  cor$frac = in_data$frac
  cor
}

run_fractional_pearson = function(in_data, log_transform = FALSE, replace_0 = FALSE){
  if (log_transform) {
    use_data = log(in_data$data)
  } else {
    use_data = in_data$data
  }
  if (replace_0) {
    tmp_data = use_data
    tmp_data[is.na(tmp_data)] = 0
    use_data = tmp_data
  }
  
  t1 = Sys.time()
  out_cor = cor(use_data, method = "pearson", use = "pairwise.complete.obs")
  t2 = Sys.time()
  t_res = as.numeric(difftime(t2, t1, units = "secs"))
  res = list(cor = out_cor,
             run_time = t_res,
             type = in_data$type,
             frac = in_data$frac)
  res
}

run_fractional_kendall = function(in_data, replace_0 = FALSE){
  if (replace_0) {
    tmp_data = in_data$data
    tmp_data[is.na(tmp_data)] = 0
    in_data$data = tmp_data
  }
  t1 = Sys.time()
  out_cor = multicore_kendallt(in_data$data)
  t2 = Sys.time()
  t_res = as.numeric(difftime(t2, t1, units = "secs"))
  res = list(cor = out_cor,
             run_time = t_res,
             type = in_data$type,
             frac = in_data$frac)
  res
}

multicore_kendallt = function(data_matrix){
  
  n_sample = ncol(data_matrix)
  # set everything to NA and let R take care of it
  
  ncore = future::nbrOfWorkers()
  names(ncore) = NULL
  
  pairwise_comparisons = combn(n_sample, 2)
  
  n_todo = ncol(pairwise_comparisons)
  n_each = ceiling(n_todo / ncore)
  
  split_comparisons = vector("list", ncore)
  start_loc = 1
  
  for (isplit in seq_along(split_comparisons)) {
    stop_loc = min(start_loc + n_each, n_todo)
    
    split_comparisons[[isplit]] = pairwise_comparisons[, start_loc:stop_loc, drop = FALSE]
    start_loc = stop_loc + 1
    
    if (start_loc > n_todo) {
      break()
    }
  }
  
  null_comparisons = purrr::map_lgl(split_comparisons, is.null)
  split_comparisons = split_comparisons[!null_comparisons]
  
  do_split = function(do_comparisons, data_matrix) {
    #seq_range = seq(in_range[1], in_range[2])
    #print(seq_range)
    tmp_cor = matrix(0, nrow = ncol(data_matrix), ncol = ncol(data_matrix))
    rownames(tmp_cor) = colnames(tmp_cor) = colnames(data_matrix)
    
    for (icol in seq(1, ncol(do_comparisons))) {
      iloc = do_comparisons[1, icol]
      jloc = do_comparisons[2, icol]
      tmp_cor[iloc, jloc] = tmp_cor[jloc, iloc] = cor(x = data_matrix[, iloc], y = data_matrix[, jloc], use= "pairwise.complete.obs", method = "kendall")
    }
    
    tmp_cor
  }
  # we record how much time is actually spent doing ICI-Kt
  # itself, as some of the other operations will add a bit of time
  # 
  t1 = Sys.time()
  split_cor = furrr::future_map(split_comparisons, do_split, data_matrix)
  t2 = Sys.time()
  t_diff = as.numeric(difftime(t2, t1, units = "secs"))
  
  cor_matrix = matrix(0, nrow = ncol(data_matrix), ncol = ncol(data_matrix))
  rownames(cor_matrix) = colnames(cor_matrix) = colnames(data_matrix)
  for (isplit in split_cor) {
    cor_matrix = cor_matrix + isplit
  }
  
  cor_matrix
}
