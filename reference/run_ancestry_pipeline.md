# lantern: High-Level Workflow Functions

Convenience functions for common ancestry-specific analysis workflows.
Run ancestry splitting pipeline

## Usage

``` r
run_ancestry_pipeline(gt_matrix, pt_matrix, verbose = TRUE)
```

## Arguments

- gt_matrix:

  Integer matrix of genotypes (rows=variants, cols=samples) Values:
  0=homozygous ref, 1=heterozygous, 2=homozygous alt rownames: variant
  IDs (e.g., "chr:pos" or any identifier) colnames: sample IDs

- pt_matrix:

  Integer matrix of parent-of-origin ancestry codes (rows=samples,
  cols=regions) Values: 1=EUR/EUR, 2=AFR/EUR (mixed), 3=AFR/AFR
  rownames: sample IDs (must match colnames of gt_matrix) colnames:
  region IDs (e.g., "chr:start-end")

- verbose:

  Print progress messages (default TRUE)

## Value

List with elements:

- african:

  African ancestry-specific dosage matrix

- european:

  European ancestry-specific dosage matrix

- counts:

  List of ancestry counts per region

- overlap:

  Overlap information

## Details

Complete pipeline: split genotype matrix by ancestry, return dosage
matrices. Handles sample and variant mismatches by using only
overlapping data. Memory efficient: single-pass subsetting using index
vectors.

## Examples

``` r
# Example with mismatched dimensions
gt <- matrix(c(2, 1, 0, 1, 2, 1), nrow = 2, ncol = 3,
             dimnames = list(c("chr1:100", "chr1:200"),
                            c("sample_A", "sample_B", "sample_C")))

pt <- matrix(c(3, 2, 1, 3, 2, 1), nrow = 3, ncol = 2,
             dimnames = list(c("sample_A", "sample_B", "sample_D"),
                            c("chr1:50-150", "chr1:150-250")))

result <- run_ancestry_pipeline(gt, pt)
#> === LANTERN Ancestry Pipeline ===
#> Step 1: Finding sample overlap...
#>   Samples: 3 in GT, 3 in PT
#>   Common samples: 2
#>   Dropped: sample_C, sample_D
#> 
#> Step 2: Finding variant/region overlap...
#>   No exact matches. Trying coordinate-based matching...
#>   Found 2 variants in 2 regions
#>   Variants in GT: 2
#>   Regions in PT: 2
#>   Common variants/regions: 4
#> 
#> Step 3: Subsetting matrices...
#> 
#> Step 3b: Filtering monomorphic variants...
#> 
#> Step 4: Splitting genotypes by ancestry...
#>   -> African dosage matrix: 2 x 2
#>   -> European dosage matrix: 2 x 2
#> 
#> Step 5: Counting ancestries per region...
#>   -> African (3): median = 1
#>   -> European (1): median = 0
#>   -> Mixed (2): median = 1
#> 
#> === Pipeline Complete ===
# Only sample_A and sample_B are kept (common to both)
# Only chr1:100 and chr1:200 that overlap regions are kept
```
