# Split unphased genotype matrix by ancestry for K populations

Generalises
[`split_by_ancestry`](https://yuw444.github.io/LANTERN/reference/split_by_ancestry.md)
to an arbitrary number of ancestries. For each variant, population
allele proportions \\p_1, \ldots, p_K\\ (summing to 1) are estimated
from samples with unambiguous ancestry (pure-ancestry homozygotes and
heterozygotes, plus mixed hom-alts which each contribute one allele to
each parent pool). Ambiguous heterozygotes in each mixed pair are then
split by the conditional pairwise ratio \\p_i / (p_i + p_j)\\. When no
unambiguous alt alleles exist for a pair (singleton edge case), the
split defaults to 0.5 / 0.5.

## Usage

``` r
split_by_ancestry_multi(gt_genotype, ancestry, pure_codes, mixed_codes)
```

## Arguments

- gt_genotype:

  Integer matrix (variants x samples) of genotypes (0/1/2).

- ancestry:

  Integer matrix (variants x samples) of ancestry codes, same dimensions
  as `gt_genotype`.

- pure_codes:

  Named integer vector mapping each population label to its
  pure-ancestry diploid code. Example:
  `c(AFR = 3L, EUR = 1L, NAT = 4L)`.

- mixed_codes:

  Data frame with three columns:

  code

  :   Integer ancestry code for this mixed pair.

  pop1

  :   Character name of the first parent population (must be a name in
      `pure_codes`).

  pop2

  :   Character name of the second parent population (must be a name in
      `pure_codes`).

  One row per mixed-ancestry diploid type. Example for 3 populations:
  `data.frame(code=c(2L,5L,6L), pop1=c("AFR","AFR","EUR"), pop2=c("EUR","NAT","NAT"))`

## Value

Named list of K numeric matrices (variants x samples), one per
population in `pure_codes`, in the same order. Each element is named
after the corresponding entry in `pure_codes`.

## See also

[`split_by_ancestry`](https://yuw444.github.io/LANTERN/reference/split_by_ancestry.md)
for the hardcoded 2-population (AFR + EUR) wrapper;
[`split_phased_multi`](https://yuw444.github.io/LANTERN/reference/split_phased_multi.md)
for the phased K-population analogue.

## Examples

``` r
# 3-population example: AFR=3, EUR=1, NAT=4; mixed codes 2/5/6
set.seed(1)
gt  <- matrix(sample(0:2, 20, replace = TRUE, prob = c(.8,.15,.05)),
              nrow = 4, ncol = 5)
# ancestry codes: 3=AFR/AFR, 1=EUR/EUR, 4=NAT/NAT,
#                 2=AFR/EUR, 5=AFR/NAT, 6=EUR/NAT
anc <- matrix(sample(c(1L,2L,3L,4L,5L,6L), 20, replace = TRUE),
              nrow = 4, ncol = 5)
pure  <- c(AFR = 3L, EUR = 1L, NAT = 4L)
mixed <- data.frame(code = c(2L, 5L, 6L),
                    pop1 = c("AFR", "AFR", "EUR"),
                    pop2 = c("EUR", "NAT", "NAT"))
out <- split_by_ancestry_multi(gt, anc, pure, mixed)
names(out)   # "AFR" "EUR" "NAT"
#> [1] "AFR" "EUR" "NAT"
```
