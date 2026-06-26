# Parse RFMix MSP file

Parse the RFMix MSP (local ancestry tract) file into ancestry matrices.
MSP format: two header lines followed by data rows. Line 1: population
codes (e.g., "#Subpopulation order/codes: AFR=0 EUR=1"). Line 2: column
headers. Data rows: ancestry calls per haplotype per tract.

## Usage

``` r
.parse_msp(msp_path, verbose = TRUE)
```

## Arguments

- msp_path:

  Path to MSP file (plain text or gzipped).

- verbose:

  Print progress messages.

## Value

List with: `pop_codes` (named integer vector), `sample_ids` (character),
`tract_df` (data.frame sorted by spos), `anc_hap0` (integer matrix:
samples × tracts), `anc_hap1` (integer matrix: samples × tracts).
