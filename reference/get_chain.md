# Extract a single chain as a results object

A results object stores every chain identically in `results$chains` and
keeps no estimates of its own. This function flattens one chain into a
standalone `snp_slice_results` object carrying that chain's
`map_allocation_matrix`, `map_dictionary_matrix`,
`final_allocation_matrix`, `final_dictionary_matrix`, and
`mcmc_samples`, which is how the estimates of a run are reached.
Diagnostic functions call it for you via their `chain` argument.

## Usage

``` r
get_chain(results, chain = results$best_chain)
```

## Arguments

- results:

  SNP-Slice results object

- chain:

  Index of the chain to extract (defaults to the best chain)

## Value

An `snp_slice_results` object for the requested chain

## Examples

``` r
result <- load_example_results()
chain1 <- get_chain(result, 1)
summary(chain1)
#> SNP-Slice Results Summary
#> ========================
#> 
#> Model: negative_binomial 
#> Data dimensions: 200 hosts x 96 SNPs
#> Data type: read_counts 
#> 
#> Results:
#> - Number of strains identified: 52 
#> - Number of hosts: 200 
#> - Multiplicity of infection (MOI):
#>   - Mean MOI: 2.61 
#>   - Median MOI: 1 
#>   - Range: 1 - 11 
#>   - Single infections: 105 ( 52.5 %)
#>   - Mixed infections: 95 ( 47.5 %)
#> 
#> Convergence:
#> - Iterations run: 1250 
#> - Samples retained (post-burn-in): 250 
#> - Gap Converged: No 
#> - Final log posterior: -74968.33 
#> - MAP log posterior: -74834.45 
#> - Final k*: 113 
#> - MAP k*: 113 
#> 
dim(chain1$map_allocation_matrix)
#> [1] 200  52
```
