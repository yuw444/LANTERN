#!/usr/bin/env Rscript
# Prepare test PLINK files with ancestry codes (01/02/03) format

library(data.table)

cat("=== Preparing PLINK files with ancestry codes ===\n")

# Load existing PT matrix
pt_file <- "test/data/1000g/subset/pt_matrix.tsv"
pt <- fread(pt_file)

cat("1. Loaded PT matrix:\n")
cat("   Samples:", nrow(pt), "\n")
cat("   Regions:", ncol(pt) - 1, "\n\n")

# Prepare sample info (.fam format)
# Family ID, Individual ID, Paternal ID, Maternal ID, Sex (1=male, 2=female), Phenotype
fam_data <- data.frame(
    famid = pt$sample_id,
    iid = pt$sample_id,
    patid = "0",
    matid = "0",
    sex = sample(c(1, 2), nrow(pt), replace = TRUE),
    pheno = "-9"
)
write.table(fam_data, "test/data/1000g/subset/ancestry.fam", 
            quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

cat("2. Created: ancestry.fam\n")

# Prepare variant info (.bim format)
# Chromosome, SNP ID, cM, Position, allele1 (REF), allele2 (ALT)
n_regions <- ncol(pt) - 1
bim_data <- data.frame(
    chr = 22,
    snp_id = paste0("REG", 1:n_regions),
    cm = 0,
    pos = 1:n_regions,
    allele1 = "A",
    allele2 = "T"
)
write.table(bim_data, "test/data/1000g/subset/ancestry.bim", 
            quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

cat("3. Created: ancestry.bim\n")

# Create .bed file (binary)
# ANCESTRY CODE -> PLINK GENOTYPE
# 01 (EUR/EUR) -> 2
# 02 (AFR/EUR) -> 1
# 03 (AFR/AFR) -> 0
ancestry_mat <- as.matrix(pt[, -1, with = FALSE])
genotype_mat <- ifelse(ancestry_mat == 1, 2,
                ifelse(ancestry_mat == 2, 1,
                ifelse(ancestry_mat == 3, 0, as.integer(NA))))

n_samples <- nrow(genotype_mat)
n_variants <- ncol(genotype_mat)
bytes_per_variant <- ceiling(n_samples / 4)

cat("4. Creating ancestry.bed (", n_variants, "variants x", n_samples, "samples)...\n")

con <- file("test/data/1000g/subset/ancestry.bed", "wb")

# Write magic bytes (PLINK binary v1 - SNP-major mode)
writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)

# Write genotypes in SNP-major mode
for (v in 1:n_variants) {
    variant_bytes <- rep(as.raw(0), bytes_per_variant)
    
    for (s in 1:n_samples) {
        byte_idx <- (s - 1) %/% 4 + 1
        bit_idx <- ((s - 1) %% 4) * 2
        
        gt <- genotype_mat[s, v]
        if (is.na(gt)) {
            gt_val <- 3L
        } else {
            gt_val <- as.integer(gt)
        }
        
        current <- as.integer(variant_bytes[byte_idx])
        new_val <- bitwOr(current, bitwShiftL(gt_val, bit_idx))
        variant_bytes[byte_idx] <- as.raw(new_val)
    }
    
    writeBin(variant_bytes, con)
}

close(con)
cat("   Created: ancestry.bed\n")

# Save ancestry code mapping
code_map <- data.frame(
    ancestry_code = c(1, 2, 3),
    meaning = c("EUR/EUR", "AFR/EUR", "AFR/AFR"),
    plink_gt = c(2, 1, 0)
)
write.table(code_map, "test/data/1000g/subset/ancestry_codes.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

cat("\n5. Created: ancestry_codes.txt (code mapping)\n")

# Verify
cat("\n=== Verification ===\n")
cat("File sizes:\n")
cat("  .fam:", file.info("test/data/1000g/subset/ancestry.fam")$size, "bytes\n")
cat("  .bim:", file.info("test/data/1000g/subset/ancestry.bim")$size, "bytes\n")
cat("  .bed:", file.info("test/data/1000g/subset/ancestry.bed")$size, "bytes\n\n")

cat("PLINK prefix: test/data/1000g/subset/ancestry\n")
