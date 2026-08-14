# Extract allocation information from SNP-Slice results

Extract allocation information from SNP-Slice results

## Usage

``` r
extract_allocations(results, chain = NULL)
```

## Arguments

- results:

  SNP-Slice results object

- chain:

  Which chain of a multi-chain run to use. `NULL` (default) uses the
  top-level estimates, i.e. the chain with the highest MAP log posterior
  (`results$best_chain`). Give an index to analyse a specific chain
  instead; see
  [`compare_chains`](https://plasmogenepi.github.io/snp.slicer/reference/compare_chains.md).

## Value

List containing allocation information
