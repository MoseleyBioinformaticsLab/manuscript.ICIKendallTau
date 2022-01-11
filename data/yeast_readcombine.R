library(dplyr)
yeast_readcombine = function(){
  yeast_files = dir(here::here("data", "barton_preprocessed"), pattern = "gbgout$", full.names = TRUE)
  
  yeast_data = purrr::map_dfr(yeast_files, yeast_read)
  
  yeast_info = yeast_data %>%
    dplyr::select(sample, rep, sample_rep) %>%
    unique()
  rownames(yeast_info) = NULL
  
  yeast_counts = yeast_data %>%
    dplyr::select(gene, count, sample_rep) %>%
    tidyr::pivot_wider(id_cols = c(gene, sample_rep),
                       names_from = sample_rep,
                       values_from = count)
  yeast_matrix = as.matrix(yeast_counts[, 2:ncol(yeast_counts)])
  rownames(yeast_matrix) = yeast_counts$gene
  list(counts = yeast_matrix, info = yeast_info)
}

yeast_read = function(in_file){
  read_data = read.table(in_file, header = FALSE, sep = "\t")
  no_feature = which(grepl("no_feature", read_data[[1]])) - 1
  read_data = read_data[1:no_feature, ]
  names(read_data) = c("gene", "count")
  
  id_parts = strsplit(basename(in_file), split = "_", fixed = TRUE)[[1]]
  rep_num = gsub("rep", "", id_parts[2])
  read_data$sample = id_parts[1]
  read_data$rep = rep_num
  read_data$sample_rep = paste0(id_parts[1], ".", rep_num)
  read_data
}

yeast_data = yeast_readcombine()
saveRDS(yeast_data, file = here::here("data", "yeast_counts_info.rds"))
