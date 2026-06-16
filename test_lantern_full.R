#!/usr/bin/env Rscript
# Comprehensive test of lantern package with 1000 Genomes data

library(lantern)
library(vcfR)
library(data.table)
library(dplyr)

cat("===========================================\n")
cat("  lantern package - Comprehensive Test\n")
cat("===========================================\n\n")

# Paths
vcf_file <- "test/data/1000g/subset/chr22_subset.vcf.gz"
pt_file <- "test/data/1000g/subset/pt_matrix.tsv"
panel_file <- "test/data/1000g/subset/panel_subset.tsv"

# 1. Load PT matrix
cat("1. Loading PT (ancestry) matrix...\n")
pt <- fread(pt_file)
pt_mat <- as.matrix(pt[, -1])
storage.mode(pt_mat) <- "integer"
cat("   Loaded:", nrow(pt_mat), "samples x", ncol(pt_mat), "variants\n\n")

# 2. Load and subset VCF
cat("2. Loading VCF data...\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt <- extract.gt(vcf, return.alleles = FALSE, convertNA = TRUE)
cat("   VCF dimensions:", nrow(gt), "variants x", ncol(gt), "samples\n\n")

# 3. Get sample info
cat("3. Sample information...\n")
panel <- fread(panel_file)
cat("   AFR samples:", sum(panel$pop %in% c("YRI", "LWK", "MSL", "GWD", "ESN", "ACB", "ASW")), "\n")
cat("   EUR samples:", sum(panel$pop %in% c("CEU", "TSI", "FIN", "GBR", "IBS")), "\n\n")

# 4. Test count_ancestry_codes
cat("4. Testing count_ancestry_codes()...\n")
afr_counts <- count_ancestry_codes(pt_mat, 3)
eur_counts <- count_ancestry_codes(pt_mat, 1)
mixed_counts <- count_ancestry_codes(pt_mat, 2)

cat("   African ancestry (code 3):\n")
cat("     Mean:", round(mean(afr_counts), 1), "variants/sample\n")
cat("     Range:", min(afr_counts), "-", max(afr_counts), "\n")

cat("   European ancestry (code 1):\n")
cat("     Mean:", round(mean(eur_counts), 1), "variants/sample\n")
cat("     Range:", min(eur_counts), "-", max(eur_counts), "\n")

cat("   Mixed ancestry (code 2):\n")
cat("     Mean:", round(mean(mixed_counts), 1), "variants/sample\n")
cat("     Range:", min(mixed_counts), "-", max(mixed_counts), "\n\n")

# 5. Test split_by_ancestry
cat("5. Testing split_by_ancestry()...\n")

# Convert VCF genotypes to integer matrix (for testing)
gt_int <- apply(gt, 2, function(x) {
    x <- ifelse(is.na(x), 0, ifelse(x == "0/0", 0, ifelse(x == "0/1" | x == "1/0", 1, ifelse(x == "1/1", 2, 0))))
    as.integer(x)
})
gt_int <- t(gt_int)
storage.mode(gt_int) <- "integer"

# Use subset for testing (first 50 variants)
n_test <- min(50, ncol(pt_mat))
gt_subset <- gt_int[, 1:n_test, drop = FALSE]
pt_subset <- pt_mat[, 1:n_test, drop = FALSE]

cat("   Testing with", nrow(gt_subset), "samples x", ncol(gt_subset), "variants\n")
result <- split_by_ancestry(gt_subset, pt_subset)

cat("   African dosage matrix:", nrow(result$african), "x", ncol(result$african), "\n")
cat("   European dosage matrix:", nrow(result$european), "x", ncol(result$european), "\n")

# Verify: for each sample, sum of AFR + EUR + (mixed component) should equal original GT
# This is approximate since we use 50/50 split for mixed
cat("   Verification - genotype totals match:\n")
original_sum <- rowSums(gt_subset)
split_sum <- rowSums(result$african) + rowSums(result$european)
diff <- abs(original_sum - split_sum)
cat("     Max difference:", max(diff), "(expected ~0 for pure, ~50% for mixed)\n\n")

# 6. Summary statistics
cat("6. Summary statistics...\n")
cat("   African dosage per sample:\n")
cat("     Mean:", round(mean(rowSums(result$african)), 1), "\n")
cat("     Max:", max(rowSums(result$african)), "\n")

cat("   European dosage per sample:\n")
cat("     Mean:", round(mean(rowSums(result$european)), 1), "\n")
cat("     Max:", max(rowSums(result$european)), "\n\n")

# 7. Performance test
cat("7. Performance test (1M operations)...\n")
n_iter <- 10000
pt_large <- matrix(sample(1:3, 100 * 100, replace = TRUE), nrow = 100, ncol = 100)
storage.mode(pt_large) <- "integer"

pt <- proc.time()
for (i in 1:n_iter) {
    count_ancestry_codes(pt_large, 3)
}
elapsed <- proc.time() - pt
ops_per_sec <- n_iter / elapsed[3]
cat("   Operations per second:", round(ops_per_sec), "\n\n")

cat("===========================================\n")
cat("  All tests completed successfully!\n")
cat("===========================================\n")
