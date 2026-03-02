# Summarize allele frequencies from per-sample counts (MCMC path). Builds per-sample frequencies, then returns posterior mean, SD, credible interval, and per-sample mean count. Sample-size invariant.

Summarize allele frequencies from per-sample counts (MCMC path). Builds
per-sample frequencies, then returns posterior mean, SD, credible
interval, and per-sample mean count. Sample-size invariant.

## Usage

``` r
summarize_allele_frequencies_mcmc(
  allele_counts_list,
  interval = 0.95,
  n_samples = length(allele_counts_list)
)
```
