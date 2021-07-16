do_yeast_matchpaper = function(yeast_counts){
  yeast_counts[yeast_counts == 0] = NA
  yeast_cor = cor(yeast_counts, use = "pairwise.complete.obs")
  yeast_cor
}

do_yeast_remove0 = function(yeast_counts){
  yeast_counts[yeast_counts == 0] = NA
  yeast_counts = log(yeast_counts)
  yeast_cor = cor(yeast_counts, use = "pairwise.complete.obs")
  yeast_cor
}
