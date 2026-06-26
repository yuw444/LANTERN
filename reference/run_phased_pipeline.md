# Run phased ancestry splitting pipeline

Complete pipeline for phased data: parse an RFMix MSP file and a phased
VCF, intersect samples, map each variant to its ancestry tract, call
[`split_phased`](https://yuw444.github.io/LANTERN/reference/split_phased.md)
via the C backend, and optionally write ancestry-specific VCFs with a DS
field that can be converted to GDS.

## Usage

``` r
run_phased_pipeline(
  vcf_path,
  msp_path,
  out_path,
  chrom = NULL,
  write_vcf = TRUE,
  verbose = TRUE
)
```

## Arguments

- vcf_path:

  Path to phased VCF/BCF file (plain or gzipped). `bcftools` must be in
  `PATH`.

- msp_path:

  Path to RFMix MSP file (plain text or gzipped TSV).

- out_path:

  Output directory for VCFs, GDS, and cache files.

- chrom:

  Chromosome to process (e.g., `"chr19"` or `"19"`). If `NULL`, all
  chromosomes present in the VCF are used.

- write_vcf:

  Logical; if `TRUE`, write ancestry-specific VCFs and convert them to
  GDS format using
  [`SeqArray::seqVCF2GDS()`](https://rdrr.io/pkg/SeqArray/man/seqVCF2GDS.html).

- verbose:

  Print step-by-step progress messages.

## Value

Invisibly, a list with elements:

- african:

  Numeric matrix (variants × samples) of African dosages.

- european:

  Numeric matrix (variants × samples) of European dosages.

- variant_info:

  data.frame with columns chrom, pos, ref, alt.

- sample_ids:

  Character vector of common sample IDs.

- tract_info:

  data.frame of ancestry tracts from the MSP file.

- overlap:

  List of intersection statistics.

- vcf_paths:

  (when `write_vcf = TRUE`) Named list of output file paths.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_phased_pipeline(
  vcf_path = "data/chr19.phased.bcf",
  msp_path = "data/chr19.msp.tsv.gz",
  out_path = "output/",
  chrom    = "chr19",
  write_vcf = TRUE
)
head(result$variant_info)
} # }
```
