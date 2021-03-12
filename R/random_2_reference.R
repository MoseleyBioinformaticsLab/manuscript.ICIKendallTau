random_2_reference = function(in_cor, ref_cor){
  tri_cor = extract_data(in_cor)
  tri_ref = extract_data(ref_cor)
  
  diff_cor = tri_cor - tri_ref
  data.frame(median = median(diff, na.rm = TRUE),
             sd = sd(diff, na.rm = TRUE),
             time = in_cor$run_time,
             type = in_cor$type,
             fraction = in_cor$frac)
}

extract_data = function(data){
  lower_tri = lower.tri(data$cor, diag = FALSE)
  data$cor[lower_tri]
}
