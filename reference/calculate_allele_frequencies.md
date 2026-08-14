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
  estimate = c("final_sample", "map", "posterior"),
  n_samples = 100,
  interval = 0.95,
  allele_sep = "|",
  additional_burnin = 0,
  chain = NULL
)
```

## Arguments

- results:

  A `snp_slice_results` object containing MCMC results

- snp_indices:

  A vector of SNP indices to treat as a single allele

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

- allele_sep:

  Separator for allele strings (default: "\|")

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

The structure depends on `estimate`.

- Point estimate (`"final_sample"` or `"map"`):

  Data frame with columns: `allele` (string representation of the
  allele, e.g. `"ref|alt|ref"` for 3 SNPs), `frequency` (proportion of
  total parasites with this allele; sums to 1), `count` (number of
  parasites with this allele), `total_parasites` (total parasites; same
  for every row).

- Posterior (`"posterior"`):

  Data frame with columns: `allele`, `frequency` (posterior mean
  proportion), `frequency_sd` (posterior SD of proportion across
  samples), `frequency_lower` and `frequency_upper` (credible interval,
  e.g. 2.5% and 97.5%), `mean_count` (mean parasites with this allele
  per MCMC sample; does not scale with `n_samples`), `n_samples`.
  Attribute `mean_total_parasites`: mean total parasites per MCMC sample
  (same for all alleles).

## Details

With `estimate = "posterior"`, counts are summarized as per-sample means
rather than sums, so `mean_count` and `mean_total_parasites` are
interpretable regardless of `n_samples`. Frequency uncertainty
(`frequency_sd`, `frequency_lower`, `frequency_upper`) is computed from
the distribution of allele frequencies across MCMC samples.

## Examples

``` r
# Load example results
result <- load_example_results()

# Calculate allele frequencies for SNPs 1, 5, and 10 (final sample)
allele_freqs <- calculate_allele_frequencies(result, c(1, 5, 10))
print(allele_freqs)
#>        allele   frequency count total_parasites
#> 8 ref|ref|ref 0.634285714   333             525
#> 2 ref|alt|alt 0.135238095    71             525
#> 3 alt|ref|alt 0.091428571    48             525
#> 4 ref|ref|alt 0.074285714    39             525
#> 7 alt|ref|ref 0.036190476    19             525
#> 6 ref|alt|ref 0.017142857     9             525
#> 1 alt|alt|alt 0.005714286     3             525
#> 5 alt|alt|ref 0.005714286     3             525

# Posterior: mean, SD, credible interval, and per-sample mean count
if (!is.null(get_chain(result)$mcmc_samples)) {
  allele_freqs_post <- calculate_allele_frequencies(result, c(1, 5, 10), estimate = "posterior", n_samples = 50)
  print(allele_freqs_post)
}
#>                  allele   frequency frequency_sd frequency_lower
#> ref|ref|ref ref|ref|ref 0.634097532  0.002736118    0.6310209453
#> ref|alt|alt ref|alt|alt 0.138035001  0.005609773    0.1246050259
#> alt|ref|alt alt|ref|alt 0.082391717  0.009436347    0.0665684411
#> ref|ref|alt ref|ref|alt 0.073155009  0.004944094    0.0665684411
#> alt|ref|ref alt|ref|ref 0.042251963  0.007466272    0.0322245146
#> ref|alt|ref ref|alt|ref 0.014464276  0.002814664    0.0095278989
#> alt|alt|alt alt|alt|alt 0.007971928  0.006295963    0.0004261364
#> alt|alt|ref alt|alt|ref 0.007632574  0.003088536    0.0019108694
#>             frequency_upper mean_count n_samples
#> ref|ref|ref      0.63884162     332.32        50
#> ref|alt|alt      0.14621117      72.34        50
#> alt|ref|alt      0.09871358      43.18        50
#> ref|ref|alt      0.08405791      38.34        50
#> alt|ref|ref      0.04961832      22.14        50
#> ref|alt|ref      0.01915281       7.58        50
#> alt|alt|alt      0.01869025       4.18        50
#> alt|alt|ref      0.01486616       4.00        50
```
