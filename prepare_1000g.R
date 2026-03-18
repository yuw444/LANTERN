#!/usr/bin/env Rscript
# Prepare 1000 Genomes chr22 data for ancestry-specific analysis

library(vcfR)
library(data.table)

cat("=== Downloading and preparing 1000 Genomes chr22 test data ===\n\n")

# Paths
vcf_file <- "test/data/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
panel_file <- "test/data/1000g/integrated_call_samples_v3.20130502.ALL.panel"
out_dir <- "test/data/1000g/subset"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read sample panel
cat("1. Reading sample panel...\n")
panel <- fread(panel_file, nrows = 2504, select = c(1:4), col.names = c("sample", "pop", "super_pop", "gender"))
cat("   Total samples:", nrow(panel), "\n")
cat("   Populations:", paste(unique(panel$pop), collapse=", "), "\n")

# Define AFR and EUR populations
afr_pops <- c("YRI", "LWK", "MSL", "GWD", "ESN", "ACB", "ASW")  # African
eur_pops <- c("CEU", "TSI", "FIN", "GBR", "IBS")   # European

# Get sample IDs for each ancestry
afr_samples <- panel[pop %in% afr_pops, sample]
eur_samples <- panel[pop %in% eur_pops, sample]

cat("\n2. Filtering samples:\n")
cat("   AFR samples:", length(afr_samples), "\n")
cat("   EUR samples:", length(eur_samples), "\n")
cat("   Total:", length(afr_samples) + length(eur_samples), "\n")

# Write keep files for bcftools
cat("\n3. Writing keep files...\n")
writeLines(afr_samples, file.path(out_dir, "keep_afr.txt"))
writeLines(eur_samples, file.path(out_dir, "keep_eur.txt"))
writeLines(c(afr_samples, eur_samples), file.path(out_dir, "keep_both.txt"))

# Save panel subset
panel_subset <- panel[sample %in% c(afr_samples, eur_samples)]
fwrite(panel_subset, file.path(out_dir, "panel_subset.tsv"), sep="\t")
cat("   Panel saved to:", file.path(out_dir, "panel_subset.tsv"), "\n")

# Use bcftools to subset VCF (if available)
bcftools_available <- system("which bcftools > /dev/null 2>&1") == 0

if (bcftools_available) {
    cat("\n4. Subsetting VCF with bcftools...\n")
    
    # Subset to AFR + EUR samples
    system2("bcftools", args = c(
        "view", "-s", paste(c(afr_samples, eur_samples), collapse=","),
        "-Oz", "-o", file.path(out_dir, "chr22_subset.vcf.gz"),
        vcf_file
    ))
    
    # Index
    system2("bcftools", args = c(
        "index", "-t",
        file.path(out_dir, "chr22_subset.vcf.gz")
    ))
    
    cat("   VCF subset created: chr22_subset.vcf.gz\n")
    
} else {
    cat("\n4. bcftools not available, using vcfR...\n")
    
    # Read VCF header only to get sample list
    vcf <- read.vcfR(vcf_file, verbose = FALSE)
    
    # Get sample names
    vcf_samples <- colnames(vcf)
    keep_idx <- vcf_samples %in% c(afr_samples, eur_samples)
    
    cat("   VCF has", length(vcf_samples), "samples\n")
    cat("   Keeping", sum(keep_idx), "samples\n")
}

# Create simulated ancestry matrix (PT matrix)
# In real data, this would come from local ancestry inference
# Here we simulate based on population labels
cat("\n5. Creating simulated PT (ancestry) matrix...\n")

# Get all samples in VCF order for the subset
all_keep_samples <- c(afr_samples, eur_samples)

# Create ancestry codes based on super_pop:
# 3 = AFR/AFR (pure African)
# 1 = EUR/EUR (pure European)  
# 2 = AFR/EUR (admixed) - simulated for some samples
sample_pop <- panel[match(all_keep_samples, sample), pop]

set.seed(42)
n_samples <- length(all_keep_samples)
n_variants <- 100  # Just first 100 variants for testing

# Sample assignment: some admixed
admixed_rate <- 0.15  # 15% admixed
is_admixed <- runif(n_samples) < admixed_rate
pt_codes <- ifelse(is_admixed, 2, ifelse(sample_pop %in% afr_pops, 3, 1))

# Create PT matrix (samples x variants)
pt_matrix <- matrix(rep(pt_codes, n_variants), nrow = n_samples, ncol = n_variants)

# Add some variation (random ancestry for testing)
for (i in 1:n_variants) {
    # Randomly flip some to admixed
    flip_idx <- runif(n_samples) < 0.1
    pt_matrix[flip_idx, i] <- sample(c(1,2,3), sum(flip_idx), replace = TRUE)
}

# Save PT matrix
pt_dt <- as.data.table(pt_matrix)
pt_dt <- cbind(data.table(sample_id = all_keep_samples), pt_dt)
fwrite(pt_dt, file.path(out_dir, "pt_matrix.tsv"), sep="\t")
cat("   PT matrix saved:", nrow(pt_dt), "samples x", ncol(pt_dt)-1, "variants\n")

cat("\n=== Done! ===\n")
cat("Output directory:", out_dir, "\n")
cat("Files created:\n")
cat("  - panel_subset.tsv (sample info)\n")
cat("  - keep_afr.txt / keep_eur.txt / keep_both.txt\n")
cat("  - pt_matrix.tsv (ancestry codes)\n")
if (bcftools_available) cat("  - chr22_subset.vcf.gz (VCF subset)\n")
