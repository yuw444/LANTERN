# Write VCF files with ancestry-specific dosages

Write separate VCF files for African and European ancestry-specific
variant dosages.

## Usage

``` r
write_vcf_with_ancestry(
  vcf_path,
  gt_matrix,
  ancestry_matrix,
  output_african,
  output_european
)
```

## Arguments

- vcf_path:

  Input VCF path (for header template)

- gt_matrix:

  Genotype matrix (rows=variants, cols=samples)

- ancestry_matrix:

  Ancestry code matrix (same dimensions)

- output_african:

  Output VCF path for African ancestry

- output_european:

  Output VCF path for European ancestry
