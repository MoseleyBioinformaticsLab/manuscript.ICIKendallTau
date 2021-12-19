random_2_reference = function(in_cor, ref_cor, which = "raw"){
  tri_cor = extract_data(in_cor, which = which)
  tri_ref = extract_data(ref_cor, which = which)
  
  diff_cor = tri_cor - tri_ref
  data.frame(median = median(diff_cor, na.rm = TRUE),
             sd = sd(diff_cor, na.rm = TRUE),
             time = in_cor$run_time,
             type = in_cor$type,
             fraction = in_cor$frac)
}

extract_data = function(data, which = "raw"){
  lower_tri = lower.tri(data[[which]], diag = FALSE)
  data[[which]][lower_tri]
}
