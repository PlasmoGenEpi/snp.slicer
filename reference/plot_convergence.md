# Plot convergence diagnostics

Plot convergence diagnostics

## Usage

``` r
plot_convergence(
  results,
  type = "logpost",
  additional_burnin = 0,
  chain = NULL
)
```

## Arguments

- results:

  SNP-Slice results object

- type:

  Type of plot ("logpost", "kstar", "n_strains")

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

Plot object (if ggplot2 is available)
