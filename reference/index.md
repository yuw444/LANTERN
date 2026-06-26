# Package index

## Unphased splitting

Functions for diploid genotype matrices where ancestry is summarised as
a per-variant diploid code (1 = EUR/EUR, 2 = AFR/EUR, 3 = AFR/AFR).

- [`split_by_ancestry()`](https://yuw444.github.io/LANTERN/reference/split_by_ancestry.md)
  : Split genotype matrix by ancestry
- [`run_ancestry_pipeline()`](https://yuw444.github.io/LANTERN/reference/run_ancestry_pipeline.md)
  : lantern: High-Level Workflow Functions
- [`count_ancestry_codes()`](https://yuw444.github.io/LANTERN/reference/count_ancestry_codes.md)
  : Count ancestry codes in matrix

## Phased splitting

Functions for phased data where per-haplotype local ancestry is
available (e.g. from RFMix). Splitting is deterministic.

- [`split_phased()`](https://yuw444.github.io/LANTERN/reference/split_phased.md)
  : Split phased haplotypes by ancestry
- [`run_phased_pipeline()`](https://yuw444.github.io/LANTERN/reference/run_phased_pipeline.md)
  : Run phased ancestry splitting pipeline

## VCF utilities

Write and subset ancestry-specific VCF files.

- [`write_vcf_with_ancestry()`](https://yuw444.github.io/LANTERN/reference/write_vcf_with_ancestry.md)
  : Write VCF files with ancestry-specific dosages
- [`create_ancestry_vcfs()`](https://yuw444.github.io/LANTERN/reference/create_ancestry_vcfs.md)
  : Create ancestry-specific VCF files
- [`subset_vcf_by_range()`](https://yuw444.github.io/LANTERN/reference/subset_vcf_by_range.md)
  : Subset VCF by genomic range

## Downstream helpers

Helpers for preparing inputs to SMMAT gene-level tests.

- [`calc_gene_ancestry_weights()`](https://yuw444.github.io/LANTERN/reference/calc_gene_ancestry_weights.md)
  : Calculate ancestry weights per gene

## Package

- [`lantern`](https://yuw444.github.io/LANTERN/reference/lantern-package.md)
  [`lantern-package`](https://yuw444.github.io/LANTERN/reference/lantern-package.md)
  : lantern: Ancestry-Specific Rare Variant Analysis
