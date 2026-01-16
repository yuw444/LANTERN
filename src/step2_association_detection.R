library(data.table)
library(dplyr)
library(GMMAT)
library(optparse)
library(SeqArray)
library(future)
library(future.apply)
library(MASS)

option_list <- list(
  make_option(
    c("--african_gds"),
    type = "character",
    help = "Path to African GDS file",
    metavar = "file"
  ),
  make_option(
    c("--european_gds"),
    type = "character",
    help = "Path to European GDS file",
    metavar = "file"
  ),
  make_option(
    c("--data_file"),
    type = "character",
    help = "Path to phenotype file (csv, tsv, or rds)",
    metavar = "file"
  ),
  make_option(
    c("--response_type"),
    type = "character",
    default = c("continuous", "binary", "count"),
    help = "Response type: [default %default]",
    metavar = "type"
  ),
  make_option(
    c("--kinship_rds"),
    type = "character",
    help = "Path to kinship RDS file",
    metavar = "file"
  ),
  make_option(
    c("--out_path"),
    type = "character",
    help = "Output directory",
    metavar = "dir"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

afr_gds <- opt$african_gds
eur_gds <- opt$european_gds
data_file <- opt$data_file
kinship_rds <- opt$kinship_rds
out_path <- opt$out_path

# basic validation
if (
  is.null(afr_gds) ||
    is.null(eur_gds) ||
    is.null(data_file) ||
    is.null(kinship_rds) ||
    is.null(out_path)
) {
  stop(
    "All options --african_gds, --european_gds, --data_file, --kinship_rds and --out_path must be provided."
  )
}

# read phenotype data (supports .rds or delimited text)
if (grepl("\\.rds$", data_file, ignore.case = TRUE)) {
  df_pheno <- readRDS(data_file)
} else {
  df_pheno <- fread(data_file)
}

if (!"id" %in% colnames(df_pheno)) {
  stop("Phenotype data must contain an 'id' column.")
}

pheno_ids <- df_pheno$id

gds <- seqOpen(afr_gds)
gds_ids <- seqGetData(gds, "sample.id")
seqClose(gds)

df_kinship <- readRDS(kinship_rds)
kin_ids <- colnames(df_kinship)

# find common IDs among phenotype, GDS and kinship
ids_common <- Reduce(intersect, list(pheno_ids, gds_ids, kin_ids))
if (length(ids_common) == 0) {
  stop("No overlapping IDs found among phenotype, GDS, and kinship samples.")
}
message(sprintf("Found %d overlapping IDs", length(ids_common)))

# require kinship matrix to have row/col names
if (is.null(rownames(df_kinship)) || is.null(colnames(df_kinship))) {
  stop("df_kinship must have row and column names corresponding to sample IDs.")
}

# subset and reorder phenotype data to match the kinship order
row_idx <- match(ids_common, df_pheno$id)
df_pheno <- df_pheno[row_idx, , drop = FALSE]

col_to_subset <- match(ids_common, kin_ids)
## Convert kinship to plain matrix to avoid S4 serialization issues
Kmat <- as.matrix(df_kinship[col_to_subset, col_to_subset, drop = FALSE])


# create output directory if needed
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

name_res <- colnames(df_pheno)[2]
name_covs <- colnames(df_pheno)[3:ncol(df_pheno)]

formula_to_fit <- as.formula(
  paste0(
    name_res,
    " ~ ",
    paste(name_covs, collapse = " + ")
  )
)

family_to_use <- switch(
  opt$response_type,
  continuous = gaussian(link = "identity"),
  binary = binomial(link = "logit"),
  count = poisson(link = "log"),
  categorical = stop("Categorical response type is not yet supported.")
)

# NULL Models
model0 <- glmmkin(
  formula_to_fit,
  data = df_pheno,
  kins = Kmat,
  id = "id",
  family = family_to_use
)

# AA only variants
out_aa <- SMMAT(
  model0,
  afr_gds,
  gene_group_file,
  MAF.range = c(0, 0.5),
  miss.cutoff = 1,
  method = "davies",
  is.dosage = TRUE,
  ncores = 1
)

out_ee <- SMMAT(
  model0,
  eur_gds,
  gene_group_file,
  MAF.range = c(0, 0.5),
  miss.cutoff = 1,
  method = "davies",
  is.dosage = TRUE,
  ncores = 1
)

out_obs <- SMMAT(
  model0,
  eur_gds,
  gene_group_file,
  MAF.range = c(0, 0.5),
  miss.cutoff = 1,
  method = "davies",
  is.dosage = FALSE,
  ncores = 1
)

saveRDS(
  list(aa = out_aa, ee = out_ee, observed = out_obs),
  file = file.path(
    out_path,
    paste0(
      "/smmat_results_",
      name_res,
      "_",
      Sys.Date(),
      ".rds"
    )
  )
)