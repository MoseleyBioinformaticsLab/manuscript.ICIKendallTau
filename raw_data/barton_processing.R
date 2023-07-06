library(dplyr)
suppressPackageStartupMessages(library(DESeq2))
yeast_readcombine = function(){
  yeast_files = dir(here::here("raw_data", "barton_preprocessed"), pattern = "gbgout$", full.names = TRUE)
  
  yeast_data = purrr::map(yeast_files, yeast_read) |>
    dplyr::bind_rows()
  
  yeast_info = yeast_data %>%
    dplyr::select(sample, rep, treatment) %>%
    unique()
  rownames(yeast_info) = NULL
  
  yeast_counts = yeast_data %>%
    dplyr::select(gene, count, sample) %>%
    tidyr::pivot_wider(id_cols = gene,
                       names_from = sample,
                       values_from = count)
  yeast_matrix = as.matrix(yeast_counts[, 2:ncol(yeast_counts)])
  rownames(yeast_matrix) = yeast_counts$gene
  yeast_matrix = yeast_matrix[, yeast_info$sample]
  
  deseq_counts = DESeqDataSetFromMatrix(yeast_matrix, yeast_info, design = ~ treatment)
  deseq_counts = estimateSizeFactors(deseq_counts)
  deseq_norm = counts(deseq_counts, normalized = TRUE)
  
  has_1 = rowSums(deseq_norm > 0) > 0
  deseq_norm = deseq_norm[has_1, ]
  
  list(counts = deseq_norm, info = yeast_info)
}

yeast_read = function(in_file){
  read_data = read.table(in_file, header = FALSE, sep = "\t")
  no_feature = which(grepl("no_feature", read_data[[1]])) - 1
  read_data = read_data[1:no_feature, ]
  names(read_data) = c("gene", "count")
  
  id_parts = strsplit(basename(in_file), split = "_", fixed = TRUE)[[1]]
  rep_num = gsub("rep", "", id_parts[2])
  read_data$treatment = id_parts[1]
  read_data$rep = rep_num
  read_data$sample = paste0(id_parts[1], ".", rep_num)
  read_data
}

yeast_data = yeast_readcombine()
saveRDS(yeast_data, file = here::here("data", "barton_yeast_counts_info.rds"))
