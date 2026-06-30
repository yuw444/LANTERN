# Split phased haplotypes into K population-specific dosage matrices

Generalises
[`split_phased`](https://yuw444.github.io/LANTERN/reference/split_phased.md)
to an arbitrary number of ancestries. Each haplotype's allele is routed
to the dosage pool whose code matches its local ancestry call.
Haplotypes with unrecognised codes contribute 0 to all pools. NA
genotypes are treated as 0 (reference allele).

## Usage

``` r
split_phased_multi(
  gt_hap0,
  gt_hap1,
  anc_hap0,
  anc_hap1,
  pop_codes = c(AFR = 0L, EUR = 1L)
)
```

## Arguments

- gt_hap0:

  Integer matrix (variants x samples) of haplotype-0 alleles (0/1).

- gt_hap1:

  Integer matrix (variants x samples) of haplotype-1 alleles (0/1).

- anc_hap0:

  Integer matrix (variants x samples) of haplotype-0 ancestry codes.

- anc_hap1:

  Integer matrix (variants x samples) of haplotype-1 ancestry codes.

- pop_codes:

  Named integer vector mapping population label to ancestry code, e.g.
  `c(AFR = 0L, EUR = 1L, NAT = 2L)`. Codes must match those used in the
  MSP file (RFMix default: 0-based integers).

## Value

Named list of K numeric matrices (variants x samples), one per entry in
`pop_codes`, in the same order. List names equal the names of
`pop_codes`.

## Details

For the two-population case `split_phased_multi` and `split_phased` are
equivalent; prefer `split_phased_multi` for new code that may later be
extended to three or more ancestries.

## See also

[`split_phased`](https://yuw444.github.io/LANTERN/reference/split_phased.md)
for the two-population convenience wrapper,
[`split_by_ancestry`](https://yuw444.github.io/LANTERN/reference/split_by_ancestry.md)
for unphased (diploid-code) splitting.

## Examples

``` r
# Three-population panel: AFR=0, EUR=1, NAT=2
gt0 <- matrix(c(1L, 0L, 0L, 1L, 1L, 0L), nrow = 3, ncol = 2)
gt1 <- matrix(c(0L, 1L, 1L, 0L, 0L, 1L), nrow = 3, ncol = 2)
a0  <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 3, ncol = 2)
a1  <- matrix(c(1L, 2L, 0L, 2L, 0L, 1L), nrow = 3, ncol = 2)
out <- split_phased_multi(gt0, gt1, a0, a1, c(AFR = 0L, EUR = 1L, NAT = 2L))
names(out)   # "AFR" "EUR" "NAT"
#> [1] "AFR" "EUR" "NAT"
out$NAT      # native-ancestry dosage matrix
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    1    0
#> [3,]    0    0
```
