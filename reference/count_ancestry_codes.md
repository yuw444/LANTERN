# Count ancestry codes in matrix

Count occurrences of a specific ancestry code in each row of a matrix.

## Usage

``` r
count_ancestry_codes(mat, code)
```

## Arguments

- mat:

  Integer matrix with ancestry codes (1, 2, 3) - Rows represent genomic
  regions/windows - Columns represent samples

- code:

  Ancestry code to count: - 1 = EUR/EUR (Pure European) - 2 = AFR/EUR
  (Mixed) - 3 = AFR/AFR (Pure African)

## Value

Integer vector with count of `code` per row

## Examples

``` r
# PT matrix: 3 regions x 4 samples
pt <- matrix(c(3, 2, 1, 3, 2, 2, 1, 1, 3, 3, 2, 1), nrow = 3, ncol = 4)
count_ancestry_codes(pt, code = 3)  # Count African ancestry
#> [1] 3 0 1
```
