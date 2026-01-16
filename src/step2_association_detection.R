library(data.table)
library(dplyr)
library(GMMAT)
library(optparse)
library(SeqArray)
library(future)
library(future.apply)
library(MASS)

df_kinship <- readRDS("/scratch/g/pauer/Yu/smmat/rawdata/kinship.RDS")
full_patient_ids <- colnames(df_kinship)

col_to_subset <- match(patient_ids, full_patient_ids)
df_kinship_subset <- df_kinship[col_to_subset, col_to_subset, drop = FALSE]

## Convert kinship to plain matrix to avoid S4 serialization issues
Kmat <- as.matrix(df_kinship_subset)

# NULL Models
model0 <- glmmkin(
  y ~ 1,
  data = df_pheno,
  kins = Kmat,
  id = "id",
  family = gaussian(link = "identity")
)

# AA only variants
out_aa <- SMMAT(
  model0,
  afr_gds_file,
  gene_group_file,
  MAF.range = c(0, 0.5),
  miss.cutoff = 1,
  method = "davies",
  is.dosage = TRUE,
  ncores = 1
)

out_ee <- SMMAT(
  model0,
  eur_gds_file,
  gene_group_file,
  MAF.range = c(0, 0.5),
  miss.cutoff = 1,
  method = "davies",
  is.dosage = TRUE,
  ncores = 1
)

out_obs <- SMMAT(
  model0,
  eur_gds_file,
  gene_group_file,
  MAF.range = c(0, 0.5),
  miss.cutoff = 1,
  method = "davies",
  is.dosage = FALSE,
  ncores = 1
)
