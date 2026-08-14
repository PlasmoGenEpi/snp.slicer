# Load Example Analysis Results

Loads pre-computed SNP-Slice analysis results from the example data.
These results were generated using the negative binomial model with 3
chains, each with 1000 burn-in iterations followed by 250 retained
samples, so that
[`convergence_diagnostics`](https://plasmogenepi.github.io/snp.slicer/reference/convergence_diagnostics.md)
can report between-chain R-hat and ESS.

## Usage

``` r
load_example_results()
```

## Value

A `snp_slice_results` object containing the analysis results

## Examples

``` r
# Load the pre-computed results
result <- load_example_results()
print(result)
#> SNP-Slice Results
#> ================
#> Model: negative_binomial 
#> Dimensions: 200 hosts x 96 strains x 96 SNPs
#> Strains identified: 52 
#> Chains: 3 (best chain: 1 )
#> Gap Converged: No 
#> 
```
