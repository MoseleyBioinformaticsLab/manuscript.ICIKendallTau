library(recount)
library(dplyr)
suppressPackageStartupMessages(library(DESeq2))
load(here::here("raw_data/recount_adenocarcinoma/rse_gene_lung.Rdata"))

gene_counts = assays(rse_gene)$counts

sample_data = colData(rse_gene)

# some of these are going to be specific to the TCGA data, I think
sample_info = data.frame(project = sample_data$project,
                         sample_id = sample_data$gdc_file_id,
                         gender = sample_data$gdc_cases.demographic.gender,
                         project_name = sample_data$gdc_cases.project.name,
                         race = sample_data$gdc_cases.demographic.race,
                         sample_type = sample_data$gdc_cases.samples.sample_type,
                         primary_site = sample_data$gdc_cases.project.primary_site,
                         tumor_stage = sample_data$gdc_cases.diagnoses.tumor_stage,
                         disease_type = sample_data$cgc_file_disease_type) |>
  dplyr::mutate(treatment = dplyr::case_when(
    grepl("Tumor", sample_type) ~ "Tumor",
    grepl("Normal", sample_type) ~ "Normal"
  ))

adeno_info = dplyr::filter(sample_info, 
                           disease_type %in% "Lung Adenocarcinoma",
                           !(tumor_stage %in% c("not reported", "stage ia", "stage i", "stage ib")),
                           !(sample_type %in% c("Recurrent Tumor"))) %>%
  dplyr::mutate(sample = toupper(sample_id))

gene_counts = gene_counts[, adeno_info$sample]

deseq_counts = DESeqDataSetFromMatrix(gene_counts, adeno_info, design = ~ treatment)
deseq_counts = estimateSizeFactors(deseq_counts)
deseq_norm = counts(deseq_counts, normalized = TRUE)
# filter so we don't have any genes with zero across all samples
has_1 = rowSums(deseq_norm > 0) > 0
deseq_norm = deseq_norm[has_1, ]

out_count_info = list(counts = deseq_norm, info = adeno_info)

saveRDS(out_count_info, file = here::here("data/recount_adenocarcinoma_count_info.rds"))
