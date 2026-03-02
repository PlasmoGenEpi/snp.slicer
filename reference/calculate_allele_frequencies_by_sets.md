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
  interval = 0.95,
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

- interval:

  Credible interval width when `use_map = FALSE` (e.g. 0.95).

- allele_sep:

  Separator for allele strings (default: "\|").

## Value

A named list of data frames, one per target set. List names come from
`names(target_sets)` or `"set_1"`, `"set_2"`, etc. Each data frame has
the same structure as the return value of
[`calculate_allele_frequencies`](https://plasmogenepi.github.io/snp.slicer/reference/calculate_allele_frequencies.md):
with MAP, columns `allele`, `frequency`, `count`, `total_parasites`;
with MCMC, columns `allele`, `frequency`, `frequency_sd`,
`frequency_lower`, `frequency_upper`, `mean_count`, `n_samples`, and
attribute `mean_total_parasites`. See that function's help for the
meaning of each column.

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
