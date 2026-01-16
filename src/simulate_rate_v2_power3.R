library(data.table)
library(dplyr)
library(GMMAT)
library(optparse)
library(SeqArray)
library(future)
library(future.apply)
library(MASS)

## -------- CLI options --------
option_list <- list(
  make_option(
    c("-i", "--input_path"),
    type = "character",
    default = "/scratch/g/pauer/Yu/smmat/src/python_split/output/",
    help = "Path to the input gds files [default %default]"
  ),
  make_option(
    c("-o", "--output_path"),
    type = "character",
    default = "/scratch/g/pauer/Yu/smmat/src/simulation4/output/",
    help = "Output directory [default %default]"
  ),
  make_option(
    c("-e", "--effect_size"),
    type = "numeric",
    default = 2,
    help = "Effect size for causal variants [default %default]"
  ),
  make_option(
    c("-g", "--sigma_g"),
    type = "numeric",
    default = 1.0,
    help = "Genetic variance component sigma_g from kinship matrix [default %default]"
  ),
  make_option(
    c("-s", "--sigma_epsilon"),
    type = "numeric",
    default = 1.0,
    help = "Residual variance sigma_epsilon [default %default]"
  ),
  make_option(
    c("-p", "--rate"),
    type = "numeric",
    default = 0.6,
    help = "Causal percentage of SNP in each GENE [default %default]"
  ),
  make_option(
    c("G", "--gene_group"),
    type = "character",
    default = "/scratch/g/pauer/Yu/smmat/src/simulation4/meta/gene_groups/GOLGA6L9.tsv",
    help = "Path to gene group file [default %default]"
  ),
  make_option(
    c("A", "--african_gds"),
    type = "character",
    default = "/scratch/g/pauer/Yu/smmat/src/simulation4/meta/gene_gds/GOLGA6L9_aa.gds",
    help = "Path to African GDS file [default %default]"
  ),
  make_option(
    c("E", "--european_gds"),
    type = "character",
    default = "/scratch/g/pauer/Yu/smmat/src/simulation4/meta/gene_gds/GOLGA6L9_ee.gds",
    help = "Path to European GDS file [default %default]"
  ),
  make_option(
    c("-c", "--causal_ancestry"),
    type = "character",
    default = "AFRAFR",
    help = "Ancestry in which causal variants are simulated from AFRAFR, EUREUR, MIXED [default %default]"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
out_path <- opt$output_path
in_path <- opt$input_path
gene_group_file <- opt$gene_group
afr_gds_file <- opt$african_gds
eur_gds_file <- opt$european_gds
causal_ancestry <- opt$causal_ancestry

sigma_g <- opt$sigma_g
sigma_epsilon <- opt$sigma_epsilon

rate <- opt$rate
effect_size <- opt$effect_size

gds_aa <- seqOpen(afr_gds_file)
df_aa <- data.frame(
  CHROM = paste0("chr", seqGetData(gds_aa, "chromosome")),
  POS = seqGetData(gds_aa, "position"),
  stringsAsFactors = FALSE
)
seqClose(gds_aa)

## -------- Build variant groups and phenotypes once in main --------
gds_aa <- seqOpen(afr_gds_file)

# dosage matrices
ds_aa <- seqGetData(gds_aa, "annotation/format/DS")

seqClose(gds_aa)

gds_ee <- seqOpen(eur_gds_file)

# dosage matrices
ds_ee <- seqGetData(gds_ee, "annotation/format/DS")

seqClose(gds_ee)

## -------- Inputs --------
# Get sample IDs from GDS files and align dosage matrices to patient_ids
gds <- seqOpen(afr_gds_file)
patient_ids <- seqGetData(gds, "sample.id")
seqClose(gds)
df_kinship <- readRDS("/scratch/g/pauer/Yu/smmat/rawdata/kinship.RDS")
full_patient_ids <- colnames(df_kinship)

col_to_subset <- match(patient_ids, full_patient_ids)
df_kinship_subset <- df_kinship[col_to_subset, col_to_subset, drop = FALSE]

## Convert kinship to plain matrix to avoid S4 serialization issues
Kmat <- as.matrix(df_kinship_subset)

n_snps <- nrow(df_aa)

if (causal_ancestry == "AFRAFR") {
  ds_aa[which(!(ds_aa %in% c(0, 1, 2)))] <- 0
} else if (causal_ancestry == "EUREUR") {
  ds_ee[which(!(ds_ee %in% c(0, 1, 2)))] <- 0
} else {
  TRUE
}


library(doParallel)
registerDoParallel(cores = 20)

seeds <- 1:1000

test_rest <- foreach(s = seeds) %dopar%
  {
    set.seed(s)
    if (
      file.exists(file.path(
        out_path,
        paste0(
          "smmat_effectsize_",
          effect_size,
          "_rate_",
          rate,
          "_seed_",
          s,
          ".rds"
        )
      ))
    ) {
      return(s)
    }
    q_vec_aa <- sample(
      c(1, 0),
      size = n_snps,
      replace = TRUE,
      prob = c(rate, 1 - rate)
    )
    q_vec_ee <- sample(
      c(1, 0),
      size = n_snps,
      replace = TRUE,
      prob = c(rate, 1 - rate)
    )
    q_vec_obs <- sample(
      c(1, 0),
      size = n_snps,
      replace = TRUE,
      prob = c(rate, 1 - rate)
    )
    sigma_effect <- effect_size / sqrt(2 / pi)
    b1_aa <- abs(rnorm(n_snps, 0, sigma_effect)) * q_vec_aa
    b1_ee <- abs(rnorm(n_snps, 0, sigma_effect)) * q_vec_ee
    b1_obs <- abs(rnorm(n_snps, 0, sigma_effect)) * q_vec_obs

    mu <- if (causal_ancestry == "AFRAFR") {
      c(ds_aa %*% b1_aa)
    } else if (causal_ancestry == "EUREUR") {
      c(ds_ee %*% b1_ee)
    } else if (causal_ancestry == "MIXED") {
      c((ds_aa + ds_ee) %*% b1_obs)
    } else {
      stop("Invalid causal ancestry option.")
    }

    # Full phenotype with different causal components
    # perf2 <- peakRAM::peakRAM(
    #   {
    df_pheno <- data.frame(
      id = patient_ids,
      y = mvrnorm(
        n = 1,
        mu = mu,
        Sigma = sigma_g * Kmat + sigma_epsilon * diag(1, nrow(Kmat))
      )
    )

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
    #   }
    # )

    saveRDS(
      list(aa = out_aa, ee = out_ee, observed = out_obs),
      file = file.path(
        out_path,
        paste0(
          "smmat_effectsize_",
          effect_size,
          "_rate_",
          rate,
          "_seed_",
          s,
          ".rds"
        )
      )
    )

    return(s)
  }

print("Done all simulations.")
