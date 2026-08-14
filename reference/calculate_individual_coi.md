# Calculate estimated individual COI with uncertainty

Returns per-host complexity of infection (COI). A point estimate
(`estimate = "final_sample"` or `"map"`) gives one COI per host with
`NA` uncertainty columns; `estimate = "posterior"` computes posterior
mean, SD, and credible interval from the MCMC samples.

## Usage

``` r
calculate_individual_coi(
  results,
  estimate = c("final_sample", "map", "posterior"),
  n_samples = 100,
  interval = 0.95,
  additional_burnin = 0,
  chain = NULL
)
```

## Arguments

- results:

  A `snp_slice_results` object.

- estimate:

  Which estimate to report. One of `"final_sample"` (the final sample of
  the chain, matching the SNP-Slice paper; the default), `"map"` (the
  maximum a posteriori state), or `"posterior"` (posterior mean, SD, and
  credible interval over retained MCMC samples). `"final_sample"` and
  `"map"` are point estimates; only `"posterior"` needs
  `store_mcmc = TRUE`.

- n_samples:

  Number of MCMC samples to use when `estimate = "posterior"` (default:
  100)

- interval:

  Numeric in (0, 1). Credible interval width when
  `estimate = "posterior"` (e.g. 0.95).

- additional_burnin:

  Number of additional stored samples to discard from the start of the
  retained chain. The retained chain is already post-burn-in, so this
  defaults to 0.

- chain:

  Which chain of a multi-chain run to use. `NULL` (default) uses the
  top-level estimates, i.e. the chain with the highest MAP log posterior
  (`results$best_chain`). Give an index to analyse a specific chain
  instead; see
  [`compare_chains`](https://plasmogenepi.github.io/snp.slicer/reference/compare_chains.md).

## Value

A data frame with one row per host: host_index, host_id, coi_estimate,
coi_sd, coi_lower, coi_upper. Uncertainty columns are NA for a point
estimate or when no MCMC samples are available.

## Examples

``` r
result <- load_example_results()
coi_final <- calculate_individual_coi(result, estimate = "final_sample")
head(coi_final)
#>   host_index    host_id coi_estimate coi_sd coi_lower coi_upper
#> 1          1 specimen_1            1     NA        NA        NA
#> 2          2 specimen_2            1     NA        NA        NA
#> 3          3 specimen_3            1     NA        NA        NA
#> 4          4 specimen_4            1     NA        NA        NA
#> 5          5 specimen_5            1     NA        NA        NA
#> 6          6 specimen_6            1     NA        NA        NA
if (!is.null(get_chain(result)$mcmc_samples)) {
  coi_post <- calculate_individual_coi(result, estimate = "posterior", n_samples = 50)
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
