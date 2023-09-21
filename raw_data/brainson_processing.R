library(dplyr)
suppressPackageStartupMessages(library(DESeq2))

brainsonrnaseq_counts = readRDS(here::here("raw_data","brainson_rnaseq", "brainson_rnaseq201901_counts.rds"))
brainsonrnaseq_info = readRDS(here::here("raw_data", "brainson_rnaseq", "brainson_rnaseq201901_info.rds"))
brainsonrnaseq_info = brainsonrnaseq_info |>
  dplyr::mutate(treatment = type)
brainsonrnaseq_counts = brainsonrnaseq_counts[, brainsonrnaseq_info$sample]
  
deseq_counts = DESeqDataSetFromMatrix(round(brainsonrnaseq_counts), brainsonrnaseq_info, design = ~ treatment)
deseq_counts = estimateSizeFactors(deseq_counts)
deseq_norm = counts(deseq_counts, normalized = TRUE)
  
has_1 = rowSums(deseq_norm > 0) > 0
deseq_norm = deseq_norm[has_1, ]
  
brainson_out = list(counts = deseq_norm, info = brainsonrnaseq_info,
                    data_id = "brainsonrnaseq_egfrgenotype")

saveRDS(brainson_out, file = here::here("data", "brainsonrnaseq_type_counts_info.rds"))

brainsonrnaseq_info2 = brainsonrnaseq_info |>
  dplyr::mutate(treatment = tumor)

brainson_out2 = list(counts = deseq_norm, info = brainsonrnaseq_info2,
                     data_id = "brainsonrnaseq_egfrgenotypetumorculture")

saveRDS(brainson_out2, file = here::here("data", "brainsonrnaseq_type_counts_info_tumor.rds"))
