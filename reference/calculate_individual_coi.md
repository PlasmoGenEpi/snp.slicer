# Calculate estimated individual COI with uncertainty

Returns per-host complexity of infection (COI). With `use_map = TRUE`
you get a point estimate (MAP); with `use_map = FALSE` and MCMC samples,
posterior mean, SD, and credible interval are computed.

## Usage

``` r
calculate_individual_coi(
  results,
  use_map = TRUE,
  n_samples = 100,
  interval = 0.95
)
```

## Arguments

- results:

  A `snp_slice_results` object.

- use_map:

  If `TRUE` (default), use MAP only; uncertainty columns are `NA`. If
  `FALSE`, use MCMC samples for mean, SD, and interval.

- n_samples:

  When `use_map = FALSE`, number of MCMC samples to use (capped at
  available post-burnin samples).

- interval:

  Numeric in (0, 1). Credible interval width when using MCMC (e.g. 0.95
  for 2.5 and 97.5 percent quantiles).

## Value

A data frame with one row per host: host_index, host_id, coi_estimate,
coi_sd, coi_lower, coi_upper. Uncertainty columns are NA when using MAP
or when no MCMC samples are available.

## Examples

``` r
result <- load_example_results()
coi_map <- calculate_individual_coi(result, use_map = TRUE)
head(coi_map)
#>   host_index    host_id coi_estimate coi_sd coi_lower coi_upper
#> 1          1 specimen_1            1     NA        NA        NA
#> 2          2 specimen_2            1     NA        NA        NA
#> 3          3 specimen_3            1     NA        NA        NA
#> 4          4 specimen_4            1     NA        NA        NA
#> 5          5 specimen_5            1     NA        NA        NA
#> 6          6 specimen_6            1     NA        NA        NA
if (!is.null(result$mcmc_samples)) {
  coi_post <- calculate_individual_coi(result, use_map = FALSE, n_samples = 50)
  head(coi_post)
}
#>   host_index    host_id coi_estimate coi_sd coi_lower coi_upper
#> 1          1 specimen_1            1      0         1         1
#> 2          2 specimen_2            1      0         1         1
#> 3          3 specimen_3            1      0         1         1
#> 4          4 specimen_4            1      0         1         1
#> 5          5 specimen_5            1      0         1         1
#> 6          6 specimen_6            1      0         1         1
```
