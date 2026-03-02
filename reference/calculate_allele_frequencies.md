# Calculate Allele Frequencies from MCMC Results

## Usage

``` r
calculate_allele_frequencies(
  results,
  snp_indices,
  use_map = TRUE,
  n_samples = 100,
  interval = 0.95,
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

- interval:

  Numeric in (0, 1). Credible interval width when using MCMC (e.g.
  0.95).

- allele_sep:

  Separator for allele strings (default: "\|")

## Value

The structure depends on `use_map`.

- MAP (`use_map = TRUE`):

  Data frame with columns: `allele` (string representation of the
  allele, e.g. `"ref|alt|ref"` for 3 SNPs), `frequency` (proportion of
  total parasites with this allele; sums to 1), `count` (number of
  parasites with this allele in the MAP allocation), `total_parasites`
  (total parasites in the MAP allocation; same for every row).

- MCMC (`use_map = FALSE`):

  Data frame with columns: `allele`, `frequency` (posterior mean
  proportion), `frequency_sd` (posterior SD of proportion across
  samples), `frequency_lower` and `frequency_upper` (credible interval,
  e.g. 2.5\\

Calculates allele frequencies for a collection of SNPs treated as a
single allele. The function takes SNP indices and calculates the
frequency of each possible allele combination across all individuals
based on their strain allocations. With `use_map = FALSE`, counts are
summarized as per-sample means rather than sums, so `mean_count` and
`mean_total_parasites` are interpretable regardless of `n_samples`.
Frequency uncertainty (`frequency_sd`, `frequency_lower`,
`frequency_upper`) is computed from the distribution of allele
frequencies across MCMC samples. \# Load example results result \<-
load_example_results()# Calculate allele frequencies for SNPs 1, 5, and
10 (MAP) allele_freqs \<- calculate_allele_frequencies(result, c(1, 5,
10)) print(allele_freqs)# With MCMC: posterior mean, SD, credible
interval, and per-sample mean count if (!is.null(result\$mcmc_samples))
allele_freqs_mcmc \<- calculate_allele_frequencies(result, c(1, 5, 10),
use_map = FALSE, n_samples = 50) print(allele_freqs_mcmc)
