# Run one or more MCMC chains

Orchestrates
[`run_chain`](https://plasmogenepi.github.io/snp.slicer/reference/run_chain.md):
draws per-chain seeds, dispatches the chains (serially or in parallel),
and collects them into a single result. The chain that reaches the
highest MAP log posterior is reported as the best chain and drives the
top-level estimates.

## Usage

``` r
run_chains(..., n_chains = 1L, n_cores = 1L, verbose = TRUE)
```

## Arguments

- ...:

  Parameters passed on to `run_chain`

- n_chains:

  Number of independent chains to run

- n_cores:

  Number of cores used to run chains simultaneously

- verbose:

  Whether to print progress. Named explicitly because this function
  reports on the chain set as a whole as well as passing it down.

## Value

List with `chains` (per-chain results) and `best_chain`.
[`create_results_object()`](https://plasmogenepi.github.io/snp.slicer/reference/create_results_object.md)
keeps those two and the shared `parameters`; nothing else is carried
through to the results object.
