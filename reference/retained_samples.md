# Remove additional burn-in samples from MCMC results

Stored MCMC samples already exclude burn-in, so this returns the full
stored chain unless `additional_burnin` asks for further trimming.

## Usage

``` r
retained_samples(results, additional_burnin = 0)
```

## Arguments

- results:

  SNP-Slice results object

- additional_burnin:

  Number of additional stored samples to discard from the start of the
  retained chain

## Value

List of MCMC samples
