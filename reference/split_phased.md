# Split phased haplotypes by ancestry

Split phased haplotypes deterministically into African and European
ancestry-specific dosage matrices. Each haplotype's allele is added to
whichever ancestry pool matches its local ancestry call; haplotypes with
unrecognised ancestry codes contribute nothing to either pool.

## Usage

``` r
split_phased(
  gt_hap0,
  gt_hap1,
  anc_hap0,
  anc_hap1,
  pop_codes = c(AFR = 0L, EUR = 1L)
)
```

## Arguments

- gt_hap0:

  Integer matrix (variants × samples) of haplotype-0 alleles (0/1).

- gt_hap1:

  Integer matrix (variants × samples) of haplotype-1 alleles (0/1).

- anc_hap0:

  Integer matrix (variants × samples) of haplotype-0 ancestry codes.

- anc_hap1:

  Integer matrix (variants × samples) of haplotype-1 ancestry codes.

- pop_codes:

  Named integer vector of length 2 giving the AFR and EUR ancestry codes
  as used in the MSP file. Defaults to `c(AFR = 0L, EUR = 1L)` (RFMix
  convention).

## Value

List with two numeric matrices:

- african:

  African ancestry-specific dosage (variants × samples)

- european:

  European ancestry-specific dosage (variants × samples)

## Examples

``` r
gt0 <- matrix(c(1L, 0L, 0L, 1L), 2, 2)
gt1 <- matrix(c(0L, 1L, 1L, 0L), 2, 2)
a0  <- matrix(c(0L, 1L, 0L, 1L), 2, 2)
a1  <- matrix(c(1L, 0L, 1L, 0L), 2, 2)
split_phased(gt0, gt1, a0, a1)
#> $african
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    1    0
#> 
#> $european
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    0    1
#> 
```
