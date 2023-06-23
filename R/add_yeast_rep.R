add_yeast_rep = function(yeast_info){
  wt_snf_count = data.frame(rep_start = c(0, 48),
                            Sample = c("WT", "SNF2"))
  yeast_info2 = dplyr::left_join(yeast_info, wt_snf_count, by = "Sample")
  yeast_info2 = yeast_info2 %>%
    dplyr::mutate(base_rep = BiolRep + rep_start)
}
