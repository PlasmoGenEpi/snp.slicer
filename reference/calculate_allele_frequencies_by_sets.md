# Calculate allele frequencies for multiple target sets

For each target set, computes allele frequencies and returns a list of
frequency tables (one per set). Each set is a vector of target indices
or target names. The estimator is chosen by the `estimate` argument
passed through `...` to
[`calculate_allele_frequencies`](https://plasmogenepi.github.io/snp.slicer/reference/calculate_allele_frequencies.md).

## Usage

``` r
calculate_allele_frequencies_by_sets(results, target_sets, ...)
```

## Arguments

- results:

  A `snp_slice_results` object containing MCMC results.

- target_sets:

  List of vectors; each element is target indices (integer) or target
  names (character) defining one set. If the list is named, those names
  are used for the returned list.

- ...:

  Arguments passed on to
  [`calculate_allele_frequencies`](https://plasmogenepi.github.io/snp.slicer/reference/calculate_allele_frequencies.md).

## Value

A named list of data frames, one per target set. List names come from
`names(target_sets)` or `"set_1"`, `"set_2"`, etc. Each data frame has
the same structure as the return value of
[`calculate_allele_frequencies`](https://plasmogenepi.github.io/snp.slicer/reference/calculate_allele_frequencies.md):
for a point estimate (`"final_sample"` or `"map"`), columns `allele`,
`frequency`, `count`, `total_parasites`; for `"posterior"`, columns
`allele`, `frequency`, `frequency_sd`, `frequency_lower`,
`frequency_upper`, `mean_count`, `n_samples`, and attribute
`mean_total_parasites`. See that function's help for the meaning of each
column.

## Examples

``` r
result <- load_example_results()
target_sets <- list(locus_a = c(1, 5), locus_b = c(10))
freqs <- calculate_allele_frequencies_by_sets(result, target_sets)
print(freqs$locus_a)
#>    allele  frequency count total_parasites
#> 4 ref|ref 0.70857143   372             525
#> 2 ref|alt 0.15238095    80             525
#> 3 alt|ref 0.12761905    67             525
#> 1 alt|alt 0.01142857     6             525
```
