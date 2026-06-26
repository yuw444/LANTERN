# Split genotype matrix by ancestry

Split a genotype matrix into African and European ancestry-specific
dosage matrices based on a parent-of-origin ancestry matrix.

## Usage

``` r
split_by_ancestry(gt_genotype, ancestry)
```

## Arguments

- gt_genotype:

  Integer matrix of genotypes (0, 1, 2) - Rows represent variants -
  Columns represent samples

- ancestry:

  Integer matrix of ancestry codes (1, 2, 3) - Same dimensions as
  gt_genotype - Values: 1=EUR/EUR, 2=AFR/EUR (mixed), 3=AFR/AFR

## Value

List with two elements:

- african:

  Matrix of African ancestry-specific dosages

- european:

  Matrix of European ancestry-specific dosages

## Examples

``` r
# Genotype matrix: 5 variants x 4 samples
gt <- matrix(c(2, 1, 0, 1, 2, 1, 0, 2, 1, 0,
               1, 1, 1, 0, 2, 1, 0, 1, 0, 1), nrow = 5, ncol = 4)

# Ancestry matrix: same dimensions
ancestry <- matrix(c(3, 2, 1, 3, 2, 1, 2, 2, 1, 1,
                     3, 1, 2, 1, 3, 2, 1, 2, 1, 3), nrow = 5, ncol = 4)

result <- split_by_ancestry(gt, ancestry)
result$african
#>      [,1] [,2] [,3] [,4]
#> [1,]    2    0    1 0.75
#> [2,]    0    0    0 0.00
#> [3,]    0    1    1 1.00
#> [4,]    1    0    0 0.00
#> [5,]    1    0    2 1.00
result$european
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    1    0 0.25
#> [2,]    1    0    1 0.00
#> [3,]    0    1    1 1.00
#> [4,]    0    1    0 0.00
#> [5,]    1    0    0 0.00
```
