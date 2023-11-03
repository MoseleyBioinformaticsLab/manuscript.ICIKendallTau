ici = function(counts_info, id, keep_num, sample_col, class_col)
{
  # counts_info = tar_read(adenocarcinoma_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  # 
  # counts_info = tar_read(nsclc_counts_info)
  # keep_num = 0.25
  # sample_col = "sample"
  # class_col = "treatment"
  n_workers = future::nbrOfWorkers()
  if (n_workers == 80) {
    future::plan(multicore(workers = 50))
  }
  counts = counts_info$counts
  info = counts_info$info
  
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  # if (is.null(rownames(counts))) {
  #   rownames(counts) = paste0("f", seq_len(nrow(counts)))
  # }
  counts_filter = t(keep_non_zero_percentage(t(counts), sample_classes = info[[filter_col]],
                                             keep_num = keep_num))
  tmp_out = ici_kendalltau(counts_filter, global_na = c(NA, 0), return_matrix = FALSE)
  future::plan(multicore)
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "ici",
       full_id = id)
}

ici_completeness = function(counts_info, id, keep_num, sample_col, class_col)
{
  # counts_info = tar_read(yeast_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  n_workers = future::nbrOfWorkers()
  if (n_workers == 80) {
    future::plan(multicore(workers = 50))
  }
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
  counts_completeness = pairwise_completeness(counts_filter, return_matrix = FALSE)
  tmp_out = ici_kendalltau(counts_filter, global_na = c(NA, 0), return_matrix = FALSE)
  future::plan(multicore)
  list(cor = tmp_out,
       completeness = counts_completeness,
       data_id = counts_info$data_id,
       method_id = "ici_completeness",
       full_id = id)
}


kt = function(counts_info, id, keep_num, sample_col, class_col)
{
  # counts_info = tar_read(yeast_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  n_workers = future::nbrOfWorkers()
  if (n_workers == 80) {
    future::plan(multicore(workers = 50))
  }
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
  tmp_out = ici_kendalltau(counts_filter, global_na = c(NA), scale_max = FALSE, diag_good = FALSE, return_matrix = FALSE)
  future::plan(multicore)
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "kt",
       full_id = id)
}


pearson_base_nozero = function(counts_info, id, keep_num, sample_col, class_col)
{
  
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
  counts_filter_na = counts_filter
  counts_filter_na[counts_filter == 0] = NA
  # this one should match the Gierlinski paper values for median correlations
  tmp_out = pearson_wrapper(t(counts_filter_na))
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "pearson_base_nozero",
       full_id = id)
}

pearson_base = function(counts_info, id, keep_num, sample_col, class_col)
{
  
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
  tmp_out = pearson_wrapper(t(counts_filter))
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "pearson_base",
       full_id = id)
}

pearson_log1p = function(counts_info, id, keep_num, sample_col, class_col)
{
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
  tmp_out = pearson_wrapper(log1p(t(counts_filter)))
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "pearson_log1p",
       full_id = id)
}

pearson_log = function(counts_info, id, keep_num, sample_col, class_col)
{
  
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
  log_counts = log(counts_filter)
  log_counts[is.infinite(log_counts)] = NA
  tmp_out = pearson_wrapper(t(log_counts))
  list(cor = tmp_out,
       data_id = counts_info$data_id,
       method_id = "pearson_log",
       full_id = id)
}

pearson_wrapper = function(data_matrix)
{
  n_sample = ncol(data_matrix)
  if ("furrr" %in% utils::installed.packages()) {
    ncore = future::nbrOfWorkers()
    names(ncore) = NULL
    split_fun = furrr::future_map
  } else {
    ncore = 1
    split_fun = purrr::map
  }
  
  # generate the array of comparisons, 2 x ...,
  # where each column is a comparison between two columns of data
  message("Figuring out comparisons to do ...")
  pairwise_comparisons = utils::combn(n_sample, 2)
  
  extra_comparisons = matrix(rep(seq(1, n_sample), each = 2), nrow = 2, ncol = n_sample, byrow = FALSE)
  pairwise_comparisons = cbind(pairwise_comparisons, extra_comparisons)
  
  # create a data.frame of the comparisons by the names of the columns instead,
  # this enables indexed comparisons and named comparisons, because
  # we can have row / column names in R
  # This is now n_comparisons x 2
  named_comparisons = data.frame(s1 = colnames(data_matrix)[pairwise_comparisons[1, ]],
                                 s2 = colnames(data_matrix)[pairwise_comparisons[2, ]])
  
  
  n_todo = nrow(named_comparisons)
  
  # in my experience, the list of actual comparisons to do does not get split
  # up well across compute cores (cores, machines, etc), and we end up waiting
  # on stuff to complete.
  # Therefore, this code actually does the splitting up of comparisons across
  # the number of cores (ncore) in a list. 
  
  message("Splitting up across compute ...")
  n_each = ceiling(n_todo / ncore)
  
  which_core = rep(seq(1, ncore), each = n_each)
  which_core = which_core[1:nrow(named_comparisons)]
  
  which_core = sample(which_core, length(which_core))
  
  named_comparisons$core = which_core
  named_comparisons$raw = Inf
  named_comparisons$pvalue = Inf
  named_comparisons$taumax = Inf
  
  
  split_comparisons = split(named_comparisons, named_comparisons$core)
  
  
  # finally, this is the function for doing each set of icikt comparisons,
  # possibly across cores.
  do_split = function(do_comparisons, data_matrix) {
    #seq_range = seq(in_range[1], in_range[2])
    #print(seq_range)
    
    raw = vector("numeric", nrow(do_comparisons))
    
    
    for (irow in seq_len(nrow(do_comparisons))) {
      iloc = do_comparisons[irow, 1]
      jloc = do_comparisons[irow, 2]
      cor_res = cor(data_matrix[, iloc], data_matrix[, jloc], use = "pairwise.complete", method = "pearson")
      raw[irow] = cor_res
    }
    do_comparisons$raw = raw
    #return(ls())
    do_comparisons
  }
  # we record how much time is actually spent doing ICI-Kt
  # itself, as some of the other operations will add a bit of time
  # 
  message("Running correlations ...")
  t1 = Sys.time()
  # note here, this takes our list of comparisons, and then calls the do_split
  # function above on each of them.
  split_cor = split_fun(split_comparisons, do_split, data_matrix)
  t2 = Sys.time()
  t_diff = as.numeric(difftime(t2, t1, units = "secs"))
  
  # put all the results back together again into one data.frame
  message("Recombining results ...")
  all_cor = purrr::list_rbind(split_cor)
  rownames(all_cor) = NULL
  all_cor
}
