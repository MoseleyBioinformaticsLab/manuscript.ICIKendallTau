test_left_censorship_datasets = function(counts_info, id)
{
  # counts_info = tar_read(yeast_counts_info)
  # counts_info = tar_read(adenocarcinoma_counts_info)
  # id = "adenocarcinoma"
  # counts_info = tar_read(egfrgenotypetumorculture_counts_info)
  counts_df = tibble::as_tibble(counts_info$counts)
  if (is.null(rownames(counts_info$counts))) {
    counts_df$feature = seq_len(nrow(counts_df))
  } else {
    counts_df$feature = rownames(counts_info$counts)
  }
  
  counts_matrix = counts_df |>
    dplyr::select(-feature) |>
    as.matrix()
  
  sample_classes = counts_info$info$treatment
  
  lc_results = test_left_censorship(t(counts_matrix), sample_classes = sample_classes)
  lc_stats = lc_results$binomial_test |> broom::tidy()
  lc_stats$id = id
  n_missing = count_missingness(t(counts_matrix), sample_classes = sample_classes)
  n_missing$id = id
  
  list(lc_results = lc_results,
       stats = lc_stats,
       missingness = n_missing)
}

count_missingness = function(counts_matrix, global_na = c(NA, Inf, 0), sample_classes)
{
  counts_matrix = t(counts_matrix)
  missing_matrix = setup_missing_matrix(counts_matrix, global_na)
  counts_matrix[missing_matrix] = NA
  split_class = split(seq_len(ncol(counts_matrix)), sample_classes)
  
  n_missing = purrr::imap(split_class, \(in_samples, id){
    tmp_matrix = counts_matrix[, in_samples, drop = FALSE]
    data.frame(n_na = sum(is.na(tmp_matrix)), n = prod(dim(tmp_matrix)), treatment = id)
  }) |>
  purrr::list_rbind()
  
  n_missing = dplyr::bind_rows(n_missing,
                               data.frame(n_na = sum(n_missing$n_na),
                                          n = sum(n_missing$n),
                                          treatment = "All"))
  n_missing = n_missing |>
    dplyr::mutate(perc_na = n_na / n * 100)
  n_missing
}

#' Test for left censorship
#' 
#' Does a binomial test to check if the most likely cause of missing values
#' is due to values being below the limit of detection, or coming from a 
#' left-censored distribution.
#' 
#' @param data_matrix matrix or data.frame of numeric data
#' @param global_na what represents zero or missing?
#' @param sample_classes which samples are in which class
#' 
#' @details
#' For each feature that is missing in a group of samples, we save as a possibility
#' to test. For each sample, we calculate the median value with any missing values
#' removed. Each feature that had a missing value, we test whether the remaining 
#' non-missing values are below the sample median for those samples where the 
#' feature is non-missing. A binomial test considers the total number of features
#' instances (minus missing values) as the number of trials, and the number of
#' of features below the sample medians as the number of successes.
#' 
#' There is a bit more detail in the vignette: `vignette("testing-for-left-censorship", package = "ICIKendallTau")`
#' 
#' @seealso [ici_kendalltau()]
#' 
#' @examples
#' # this example has 80% missing due to left-censorship
#' data(missing_dataset)
#' missingness = test_left_censorship(missing_dataset)
#' missingness$values
#' missingness$binomial_test
#' 
#' @importFrom stats median
#' 
#' @export
#' @return data.frame of trials / successes, and binom.test result
test_left_censorship = function(data_matrix, 
                                global_na = c(NA, Inf, 0),
                                sample_classes = NULL)
{
  # we are using columns as samples, but manuscript has rows as samples.
  data_matrix = t(data_matrix)
  if (inherits(data_matrix, "data.frame")) {
    data_matrix = as.matrix(data_matrix)
    }
    if (is.null(sample_classes)) {
    sample_classes = rep("A", ncol(data_matrix))
  }

  split_indices = split(seq_len(ncol(data_matrix)), sample_classes)
  missing_loc = setup_missing_matrix(data_matrix, global_na)
  data_matrix_missing = data_matrix
  data_matrix_missing[missing_loc] = NA

  # split the dataset by group
  split_counts = purrr::imap(split_indices, \(in_split, split_id){
  # in_split = split_indices[[1]]

  # grab the group we want to work with
  split_missing = data_matrix_missing[, in_split, drop = FALSE]

  # count the number of missing samples for each feature,
  # and keep those that have at least one
  n_miss = rowSums(is.na(split_missing))

  if (sum(n_miss) == 0) {
  out_res = data.frame(trials = 0, success = 0, class = split_id)
  return(out_res)
  }

  keep_miss = split_missing[n_miss > 0, , drop = FALSE]

  # get sample medians
  sample_medians = calculate_matrix_medians(split_missing, use = "col", na.rm = TRUE)

  # turn the medians into a matrix to make life easier
  median_matrix = matrix(sample_medians, nrow = nrow(keep_miss),
  ncol = ncol(keep_miss), byrow = TRUE)
  # do the comparison
  keep_miss_updown = keep_miss < median_matrix

  # count how many trials we ran, and how many successes we have
  all_trials = (nrow(keep_miss_updown) * ncol(keep_miss_updown)) - sum(is.na(keep_miss_updown))
  all_success = sum(keep_miss_updown, na.rm = TRUE)

  out_res = data.frame(trials = all_trials, success = all_success, class = split_id)
  return(out_res)
  }) |>
  purrr::list_rbind()


  total_trials = sum(split_counts$trials)
  total_success = sum(split_counts$success)

  binom_res = stats::binom.test(total_success, total_trials, p = 0.5, alternative = "greater")

  if (is.null(sample_classes)) {
  total_success$class = NULL
  }

  return(list(values = split_counts,
  binomial_test = binom_res))
}


#' Calculate matrix medians
#' 
#' Given a matrix of data, calculates the median value in each column or row.
#' 
#' @param in_matrix numeric matrix of values
#' @param use character of "col" or "row" defining columns or rows
#' @param ... extra parameters to the median function
#' 
#' @export
#' @return numeric
calculate_matrix_medians = function(in_matrix, use = "col", ...)
{
  if (use %in% "row") {
  in_matrix = t(in_matrix)
  } 
  out_medians = purrr::map_dbl(seq_len(ncol(in_matrix)), \(in_col){
  stats::median(in_matrix[, in_col], ...)
  })
  return(out_medians)
}

rank_order_datasets = function(counts_info, id)
{
  # counts_info = tar_read(yeast_counts_info)
  # counts_info = tar_read(nsclc_counts_info)
  # counts_info = tar_read(egfrgenotypetumorculture_counts_info)
  counts_df = tibble::as_tibble(counts_info$counts)
  if (is.null(rownames(counts_info$counts))) {
    counts_df$feature = seq_len(nrow(counts_df))
  } else {
    counts_df$feature = rownames(counts_info$counts)
  }

  counts_matrix = counts_df |>
    dplyr::select(-feature) |>
    as.matrix()

  sample_classes = counts_info$info$treatment

  ranked_values = rank_order_data(t(counts_matrix), sample_classes = sample_classes)
  just_ranks = purrr::map(ranked_values, \(in_rank){
    in_rank$n_na_rank
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(treatment = split)

  na_values = purrr::map(ranked_values, \(in_rank){
    in_rank$ordered
  })
  n_na = purrr::imap(ranked_values, \(in_rank, id){
    use_vals = in_rank$original
    n_val = nrow(use_vals) * ncol(use_vals)
    n_na = sum(is.na(use_vals))

    data.frame(treatment = id,
                n_values = n_val,
                n_na = n_na,
                perc_na = n_na / n_val * 100)
  }) |>
    purrr::list_rbind()

  list(ranks = just_ranks,
        na_values = na_values,
        n_na = n_na,
      id = id)
}  

rank_order_data = function(data_matrix, global_na = c(NA, Inf, 0), 
                           sample_classes = NULL)
{
  data_matrix = t(data_matrix)
  if (inherits(data_matrix, "data.frame")) {
    data_matrix = as.matrix(data_matrix)
  }
  missing_loc = setup_missing_matrix(data_matrix, global_na)
  data_matrix_na = data_matrix
  data_matrix_na[missing_loc] = NA
  
  if (is.null(sample_classes)) {
    sample_classes = rep("rmf_abcd", ncol(data_matrix_na))
  }
  
  split_classes = split(colnames(data_matrix_na), sample_classes)
  
  get_ranks = function(in_na)
  {
    sample_ranks = purrr::map(seq_len(ncol(in_na)), \(in_col){
      rank(in_na[, in_col], na.last = FALSE)
    })
    sample_ranks = do.call(cbind, sample_ranks)
    
    median_rank = apply(sample_ranks, 1, median)
    n_na = rowSums(is.na(in_na))
    rank_order = order(median_rank, decreasing = TRUE)
    
    perc_missing = colSums(is.na(in_na)) / nrow(in_na)
    perc_order = order(perc_missing, decreasing = TRUE)
    list(original = in_na,
         ordered = in_na[rank_order, perc_order],
         n_na_rank = data.frame(n_na = n_na,
                                median_rank = median_rank))
  }
  
  split_ranks = purrr::imap(split_classes, \(in_split, split_id){
    split_na = data_matrix_na[, in_split]
    n_na = rowSums(is.na(split_na))
    keep_na = !(n_na == ncol(split_na))
    split_na = split_na[keep_na, ]
    
    na_info = get_ranks(split_na)
    if (!(split_id %in% "rmf_abcd")) {
      na_info$n_na_rank$split = split_id
    }
    na_info
  })
  
  if (length(split_ranks) == 1) {
    split_ranks = split_ranks[[1]]
  }
  
  return(split_ranks)
}

setup_missing_matrix = function(data_matrix, global_na)
{
  exclude_loc = matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  rownames(exclude_loc) = rownames(data_matrix)
  colnames(exclude_loc) = colnames(data_matrix)
  if (length(global_na) > 0) {
    if (any(is.na(global_na))) {
      exclude_loc[is.na(data_matrix)] = TRUE
      global_na = global_na[!is.na(global_na)]
    }
    if (any(is.infinite(global_na))) {
      exclude_loc[is.infinite(data_matrix)] = TRUE
      global_na = global_na[!is.infinite(global_na)]
    }
  }
  if (length(global_na) > 0) {
    for (ival in global_na) {
      exclude_loc[data_matrix == ival] = TRUE
    }
  }
  
  return(exclude_loc)
}
