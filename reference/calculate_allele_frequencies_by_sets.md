# Calculate allele frequencies for multiple target sets

For each target set, computes allele frequencies from MCMC/MAP results
and returns a list of frequency tables (one per set). Each set is a
vector of target indices or target names.

## Usage

``` r
calculate_allele_frequencies_by_sets(
  results,
  target_sets,
  use_map = TRUE,
  n_samples = 100,
  allele_sep = "|"
)
```

## Arguments

- results:

  A `snp_slice_results` object containing MCMC results.

- target_sets:

  List of vectors; each element is target indices (integer) or target
  names (character) defining one set. If the list is named, those names
  are used for the returned list.

- use_map:

  Logical; use MAP estimates (TRUE) or sample from MCMC (FALSE).

- n_samples:

  Number of MCMC samples to use if `use_map = FALSE` (default: 100).

- allele_sep:

  Separator for allele strings (default: "\|").

## Value

A named list of data frames. Each data frame has columns `allele`,
`frequency`, `count`, `total_parasites`. Names are from
`names(target_sets)` or `"set_1"`, `"set_2"`, etc.

## Examples

``` r
result <- load_example_results()
target_sets <- list(locus_a = c(1, 5), locus_b = c(10))
freqs <- calculate_allele_frequencies_by_sets(result, target_sets)
print(freqs$locus_a)
#>    allele  frequency count total_parasites
#> 4 ref|ref 0.70793037   366             517
#> 2 ref|alt 0.15667311    81             517
#> 3 alt|ref 0.11025145    57             517
#> 1 alt|alt 0.02514507    13             517
```
