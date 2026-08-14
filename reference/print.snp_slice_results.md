# Print SNP-Slice results

Print SNP-Slice results

## Usage

``` r
# S3 method for class 'snp_slice_results'
print(x, chain = NULL, ...)
```

## Arguments

- x:

  SNP-Slice results object

- chain:

  Which chain of a multi-chain run to use. `NULL` (default) uses the
  top-level estimates, i.e. the chain with the highest MAP log posterior
  (`results$best_chain`). Give an index to analyse a specific chain
  instead; see
  [`compare_chains`](https://plasmogenepi.github.io/snp.slicer/reference/compare_chains.md).

- ...:

  Additional arguments

## Value

Print information
