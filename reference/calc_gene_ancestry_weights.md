# Calculate ancestry weights per gene

Calculate median ancestry counts for samples within each gene's region.

## Usage

``` r
calc_gene_ancestry_weights(pt_matrix, sample_ids, gene_data)
```

## Arguments

- pt_matrix:

  PT (ancestry) matrix (rows=samples, cols=regions)

- sample_ids:

  Character vector of sample IDs

- gene_data:

  Data frame with columns: gene, chr, start, end

## Value

Data frame with gene_id, n_afr, n_eur, n_mixed columns
