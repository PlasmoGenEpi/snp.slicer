# Load Example Analysis Results

Loads pre-computed SNP-Slice analysis results from the example data.
These results were generated using the negative binomial model with 2000
MCMC iterations.

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
#> Strains identified: 49 
#> Gap Converged: No 
#> 
```
