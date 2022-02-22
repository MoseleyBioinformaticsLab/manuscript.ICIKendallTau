set_nice_widths = function(ft_table, sample_width = 0.75, outlier_width = 0.8, cor_width = 1){
  ft_width = dim(ft_table)$widths
  is_outlier = grepl(".outlier$", names(ft_width))
  is_correlation = grepl(".cor$", names(ft_width))
  is_sample = grepl("^sample", names(ft_width))
  ft_bigger = ft_table
  ft_bigger = flextable::width(ft_bigger, j = which(is_outlier), width = outlier_width)
  ft_bigger = flextable::width(ft_bigger, j = which(is_correlation), width = cor_width)
  ft_bigger = flextable::width(ft_bigger, j = which(is_sample), width = sample_width)
  ft_bigger
}

put_outlier_right = function(ft_table, outlier = ".outlier$"){
  ft_width = dim(ft_table)$widths
  is_outlier = grepl(outlier, names(ft_width))
  ft_modded = flextable::align(ft_table, j = which(is_outlier), align = "right", part = "body")
  ft_modded
}
