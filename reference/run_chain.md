# Run a single MCMC chain for SNP-Slice

This function runs a single MCMC chain for a configured snp.slicer model
object.

## Usage

``` r
run_chain(
  model_obj,
  n_sample,
  n_burnin,
  gap,
  verbose,
  store_mcmc,
  chain_id = 1L,
  seed = NULL
)
```

## Arguments

- model_obj:

  Model object

- n_sample:

  Number of post-burn-in iterations to retain

- n_burnin:

  Number of iterations discarded before sampling begins

- gap:

  Early stopping threshold

- verbose:

  Whether to print progress

- store_mcmc:

  Whether to store full MCMC samples

- chain_id:

  Integer identifier of this chain (1 when running one chain)

- seed:

  Integer seed for this chain

## Value

MCMC results

## Details

Sampling iterations are returned when `store_mcmc = TRUE`.
