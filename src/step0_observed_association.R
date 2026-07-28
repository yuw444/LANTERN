library(optparse)
option_list <- list(
  make_option(
    c("-i", "--vcf_path"),
    type = "character",
    default = NULL,
    help = "Path to phased VCF/BCF (bgzipped) with the original (unsplit) genotypes to test directly -- same file passed to step1's --vcf_path",
    metavar = "file"
  ),
  make_option(
    c("-c", "--chr_id"),
    type = "character",
    default = NULL,
    help = "Chromosome to process (must match --vcf_path's CHROM column, e.g. '22')"
  ),
  make_option(
    c("-o", "--out_path"),
    type = "character",
    default = NULL,
    help = "Output directory (created if missing); OBSERVED.gds is written here",
    metavar = "file"
  ),
  make_option(
    c("--data_file"),
    type = "character",
    default = NULL,
    help = "Path to phenotype file (csv, tsv, or rds). Must contain an 'id' column; column 2 is the response, remaining columns are covariates.",
    metavar = "file"
  ),
  make_option(
    c("--gene_group_file"),
    type = "character",
    default = NULL,
    help = "Path to gene group file (no header: gene, chr, pos, ref, alt, weight)",
    metavar = "file"
  ),
  make_option(
    c("--response_type"),
    type = "character",
    default = "continuous",
    help = "Response type: continuous, binary, or count [default: %default]",
    metavar = "type"
  ),
  make_option(
    c("--kinship_rds"),
    type = "character",
    default = NULL,
    help = "Path to kinship RDS file (square matrix, row/col names = sample IDs)",
    metavar = "file"
  ),
  make_option(
    c("--out_file"),
    type = "character",
    default = NULL,
    help = "Output TSV path for the per-gene results table (gene, p_OBSERVED columns). A sibling RDS with the same basename (extension swapped to .rds) is also written alongside it, containing list(results, smmat_results) -- smmat_results being the raw GMMAT::SMMAT() output.",
    metavar = "file"
  ),
  make_option(
    c("--ncores"),
    type = "integer",
    default = NA_integer_,
    help = "Cores for GMMAT::SMMAT(). Defaults to $SLURM_CPUS_PER_TASK when running under SLURM, else 1.",
    metavar = "n"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("vcf_path", "chr_id", "out_path", "data_file", "gene_group_file",
              "kinship_rds", "out_file")
missing_opts <- required[vapply(opt[required], is.null, logical(1))]
if (length(missing_opts) > 0) {
  stop(
    "Missing required option(s): --", paste(missing_opts, collapse = ", --"),
    "\nSee: Rscript src/step0_observed_association.R --help"
  )
}

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("GMMAT", quietly = TRUE))
  stop("Package 'GMMAT' is required (pixi run install-cran-gmmat).")
if (!requireNamespace("SeqArray", quietly = TRUE))
  stop("Package 'SeqArray' is required. Install with: BiocManager::install('SeqArray')")

dir.create(opt$out_path, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Step 1.1a: Convert the original (unsplit) VCF/BCF -- restricted to
# --chr_id -- directly to GDS via bcftools + SeqArray::seqVCF2GDS(). True
# hard-call genotypes read straight from the input file, no ancestry-split
# dosage machinery involved.
# ============================================================================
message("Converting observed (unsplit) genotypes for chr", opt$chr_id,
        " directly to GDS...")
observed_vcf_tmp <- tempfile(fileext = ".vcf")
# Filter on the CHROM field via -i rather than -t/-r: -t/-r require the
# chromosome to be declared in a ##contig header line (silently produces
# zero variants otherwise, as seen with headerless/synthetic VCFs), while
# -i matches CHROM values directly regardless of header declarations.
# Tries both with and without a leading "chr" since --chr_id's exact
# convention isn't known up front.
chrom_clean <- sub("^chr", "", opt$chr_id)
chrom_expr  <- paste0(
  "CHROM=\"", chrom_clean, "\" || CHROM=\"chr", chrom_clean, "\""
)
cmd_observed <- paste0(
  "bcftools view -i '", chrom_expr, "' ", shQuote(opt$vcf_path),
  " -Ov -o ", shQuote(observed_vcf_tmp)
)
status <- system(cmd_observed)
if (status != 0 || !file.exists(observed_vcf_tmp))
  stop("bcftools view failed while extracting chr", opt$chr_id,
       " (command: ", cmd_observed, ")")
n_observed_variants <- sum(!grepl("^#", readLines(observed_vcf_tmp)))
if (n_observed_variants == 0)
  stop("bcftools view -i '", chrom_expr, "' produced 0 variants; check that ",
       "--chr_id matches --vcf_path's CHROM column (command: ", cmd_observed, ")")

observed_gds <- file.path(opt$out_path, "OBSERVED.gds")
SeqArray::seqVCF2GDS(observed_vcf_tmp, observed_gds, verbose = FALSE)
unlink(observed_vcf_tmp)
message("Wrote: ", observed_gds,
        " (true hard-call genotypes converted directly from --vcf_path)")

# read phenotype data (supports .rds or delimited text)
if (grepl("\\.rds$", opt$data_file, ignore.case = TRUE)) {
  df_pheno <- readRDS(opt$data_file)
} else {
  df_pheno <- as.data.frame(fread(opt$data_file))
}
if (!"id" %in% colnames(df_pheno)) {
  stop("Phenotype data must contain an 'id' column.")
}

name_res <- colnames(df_pheno)[2]
name_covs <- colnames(df_pheno)[3:ncol(df_pheno)]
formula_to_fit <- as.formula(paste0(name_res, " ~ ", paste(name_covs, collapse = " + ")))

family_to_use <- switch(
  opt$response_type,
  continuous = gaussian(link = "identity"),
  binary = binomial(link = "logit"),
  count = poisson(link = "log"),
  stop("Unsupported --response_type: ", opt$response_type,
       " (use continuous, binary, or count)")
)

kinship <- readRDS(opt$kinship_rds)

ncores <- opt$ncores
if (is.na(ncores)) {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", "")
  ncores <- if (nzchar(slurm_cpus)) as.integer(slurm_cpus) else 1L
}
message("Using ncores = ", ncores, " for GMMAT::SMMAT()")

# ============================================================================
# Step 1.1b: Intersect sample IDs across pheno, kinship, and the observed
# GDS, then fit the null model and run a plain GMMAT::SMMAT()
# (is.dosage = FALSE).
# ============================================================================
gds_obs <- SeqArray::seqOpen(observed_gds)
observed_ids <- SeqArray::seqGetData(gds_obs, "sample.id")
SeqArray::seqClose(gds_obs)

ids_common <- Reduce(intersect, list(df_pheno$id, rownames(kinship), observed_ids))
if (length(ids_common) == 0)
  stop("No overlapping sample IDs found among phenotype, kinship, and OBSERVED.gds")
message("Common samples: ", length(ids_common))

pheno_sub <- df_pheno[match(ids_common, df_pheno$id), , drop = FALSE]
kin_idx   <- match(ids_common, rownames(kinship))
kin_sub   <- as.matrix(kinship[kin_idx, kin_idx, drop = FALSE])

message("Fitting null model (glmmkin)...")
model0 <- GMMAT::glmmkin(formula_to_fit, data = pheno_sub, kins = kin_sub,
                          id = "id", family = family_to_use)

message("Running GMMAT::SMMAT() on observed (unsplit) genotypes...")
smmat_observed <- GMMAT::SMMAT(model0, observed_gds, opt$gene_group_file,
                                MAF.range = c(0, 0.5), miss.cutoff = 1,
                                method = "davies", is.dosage = FALSE,
                                ncores = ncores)

results <- data.frame(
  gene       = as.character(smmat_observed$group),
  p_OBSERVED = smmat_observed$E.pval,
  stringsAsFactors = FALSE
)

fwrite(results, opt$out_file, sep = "\t")
message("Wrote: ", opt$out_file, " (per-gene p_OBSERVED table, tab-separated)")

rds_file <- paste0(tools::file_path_sans_ext(opt$out_file), ".rds")
saveRDS(list(results = results, smmat_results = smmat_observed), file = rds_file)
message(
  "Wrote: ", rds_file,
  " (list(results = <same table as the TSV above>,",
  " smmat_results = <raw GMMAT::SMMAT() output>))"
)
