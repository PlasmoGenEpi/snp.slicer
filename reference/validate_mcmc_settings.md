# Validate MCMC settings

Validate MCMC settings

## Usage

``` r
validate_mcmc_settings(n_sample, n_burnin, gap, n_chains = 1, n_cores = 1)
```

## Arguments

- n_sample:

  Number of post-burn-in iterations to retain

- n_burnin:

  Number of iterations discarded before sampling begins

- gap:

  Early stopping threshold

## Value

TRUE if valid, otherwise throws error
