# Calculate Allele Frequencies from MCMC Results

Calculates allele frequencies for a collection of SNPs treated as a
single allele. The function takes SNP indices and calculates the
frequency of each possible allele combination across all individuals
based on their strain allocations.

## Usage

``` r
calculate_allele_frequencies(
  results,
  snp_indices,
  use_map = TRUE,
  n_samples = 100,
  allele_sep = "|"
)
```

## Arguments

- results:

  A `snp_slice_results` object containing MCMC results

- snp_indices:

  A vector of SNP indices to treat as a single allele

- use_map:

  Logical, whether to use MAP estimates (TRUE) or sample from MCMC
  (FALSE)

- n_samples:

  Number of MCMC samples to use if use_map = FALSE (default: 100)

- allele_sep:

  Separator for allele strings (default: "\|")

## Value

A data frame with columns:

- allele:

  String representation of the allele (e.g., "A\|T\|T" for 3 SNPs)

- frequency:

  Proportion of total parasites with this allele

- count:

  Absolute count of parasites with this allele

- total_parasites:

  Total number of parasites across all individuals

## Examples

``` r
# Load example results
result <- load_example_results()

# Calculate allele frequencies for SNPs 1, 5, and 10
allele_freqs <- calculate_allele_frequencies(result, c(1, 5, 10))
print(allele_freqs)
#>        allele   frequency count total_parasites
#> 8 ref|ref|ref 0.642166344   332             517
#> 2 ref|alt|alt 0.148936170    77             517
#> 3 alt|ref|alt 0.073500967    38             517
#> 4 ref|ref|alt 0.065764023    34             517
#> 7 alt|ref|ref 0.036750484    19             517
#> 5 alt|alt|ref 0.025145068    13             517
#> 6 ref|alt|ref 0.007736944     4             517
#> 1 alt|alt|alt 0.000000000     0             517
```
