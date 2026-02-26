# Calculate effective sample size for MCMC diagnostics

Calculate effective sample size for MCMC diagnostics

## Usage

``` r
effective_sample_size(
  results,
  parameter = "logpost",
  method = "autocorrelation"
)
```

## Arguments

- results:

  SNP-Slice results object

- parameter:

  Parameter to analyze. Options include:

  - "logpost": Log posterior probability

  - "kstar": Number of active strains

  - "n_strains": Number of strains with non-zero allocations

  - "ktrunc": Truncation level

  - "mu": Stick-breaking weights (returns ESS for each weight)

  - "A": Allocation matrix (returns overall ESS)

  - "D": Dictionary matrix (returns overall ESS)

  - "all": Calculate ESS for all available parameters

- method:

  Method for ESS calculation. Options:

  - "autocorrelation": Standard autocorrelation-based ESS (default)

  - "batch_means": Batch means method

  - "spectral": Spectral density method

## Value

Effective sample size(s) and diagnostic information
