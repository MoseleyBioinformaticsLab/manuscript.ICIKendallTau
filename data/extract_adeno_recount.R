# downloaded http://duffel.rail.bio/recount/v2/TCGA/rse_gene_lung.Rdata

library(recount)
library(dplyr)
load(here::here("data/rse_gene_lung.Rdata"))

rse_scaled = scale_counts(rse_gene)
# we use scaled counts because we want PCA and variance calculations to "just work"
scaled_counts = assays(rse_scaled)$counts
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
                         disease_type = sample_data$cgc_file_disease_type)

adeno_info = dplyr::filter(sample_info, 
                           disease_type %in% "Lung Adenocarcinoma",
                           !(tumor_stage %in% c("not reported", "stage ia", "stage i", "stage ib")),
                           !(sample_type %in% c("Recurrent Tumor"))) %>%
  dplyr::mutate(sample_id2 = toupper(sample_id))

adeno_counts = scaled_counts[, adeno_info$sample_id2]

# filter so we don't have any genes with zero across all samples
has_1 = rowSums(adeno_counts > 0) > 0
adeno_counts = adeno_counts[has_1, ]

saveRDS(adeno_counts, file = here::here("data/recount_adeno_counts.rds"))
saveRDS(adeno_info, file = here::here("data/recount_adeno_info.rds"))
