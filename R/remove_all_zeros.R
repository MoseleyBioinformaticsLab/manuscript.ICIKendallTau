remove_all_zeros = function(count_matrix){
  n_cols = ncol(count_matrix)
  all_zero = apply(count_matrix, 1, function(.x){
    sum(.x == 0) == n_cols
  })
  count_matrix = count_matrix[!all_zero, ]
  count_matrix
}
