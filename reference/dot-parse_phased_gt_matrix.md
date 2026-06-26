# Parse phased GT matrix into haplotype integer matrices (vectorised)

Parse phased GT matrix into haplotype integer matrices (vectorised)

## Usage

``` r
.parse_phased_gt_matrix(
  gt_mat,
  n_variants,
  n_samples,
  sample_order,
  verbose = FALSE
)
```

## Arguments

- gt_mat:

  Character matrix (variants × samples).

- n_variants:

  Number of variants.

- n_samples:

  Number of samples.

- sample_order:

  Column names for output matrices.

- verbose:

  Print progress.

## Value

List: `hap0`, `hap1` (integer matrices), `missing` (logical matrix).
