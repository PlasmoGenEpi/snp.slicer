# Bayesian Nonparametric Resolution of Multi-Strain Infections

SNP-Slice is a Bayesian nonparametric method for resolving multi-strain
infections using slice sampling with stick-breaking construction. The
algorithm simultaneously unveils strain haplotypes and links them to
hosts from sequencing data.

## Usage

``` r
snp_slice(
  data,
  model = "negative_binomial",
  n_sample = 10000,
  n_burnin = NULL,
  alpha = 2.6,
  rho = if (model == "categorical") 0.5 else NULL,
  threshold = 0.001,
  gap = NULL,
  n_chains = 1,
  n_cores = 1,
  seed = NULL,
  verbose = TRUE,
  log_performance = FALSE,
  store_mcmc = FALSE,
  ...
)

snp_slice_categorical(data, e1 = 0.05, e2 = 0.05, ...)

snp_slice_poisson(data, ...)

snp_slice_binomial(data, ...)

snp_slice_negative_binomial(data, ...)

snp_slice_multinomial(data, dict_prior = "empirical", ...)
```

## Arguments

- data:

  Input data. Can be a matrix, data.frame, or file path. For read count
  data, should be a list with elements `read1` and `read0` (or `total`).
  For categorical data, can be a matrix with values 0, 0.5, or 1; or a
  long-format data.frame with columns `specimen_id`, `target_id`,
  `target_value`, and `target_count`. For a categorical data.frame,
  counts are converted to categories: ref-only -\> 0, alt-only -\> 1,
  both present -\> 0.5, zero total -\> NA. Matrix and categorical file
  inputs (e.g. `*_cat.txt`) remain supported.

- model:

  Observation model to use. One of `"categorical"`, `"poisson"`,
  `"binomial"`, `"negative_binomial"` (default), or `"multinomial"`.
  Only the multinomial model keeps targets with more than two alleles;
  the others drop them. See ‘The multinomial model’ below.

- n_sample:

  Number of post-burn-in iterations to retain (default: 10000). Burn-in
  iterations are additional: the chain runs `n_burnin + n_sample`
  iterations in total and only the last `n_sample` are retained.

- n_burnin:

  Number of iterations discarded before sampling begins. If NULL,
  defaults to `floor(n_sample / 2)`.

- alpha:

  IBP concentration parameter (default: 2.6).

- rho:

  Dictionary sparsity parameter (default: 0.5 for categorical model and
  NULL otherwise, which means use the global minor allele frequency)

- threshold:

  Threshold for identifying single infections (default: 0.001).

- gap:

  Early stopping threshold. If NULL, runs for all `n_burnin + n_sample`
  iterations.

- n_chains:

  Number of independent MCMC chains to run (default: 1). Each chain is
  seeded separately; the chain reaching the highest MAP log posterior
  supplies the top-level estimates and all chains are kept in
  `result$chains`.

- n_cores:

  Number of cores used to run chains simultaneously (default: 1, i.e.
  chains run sequentially). Capped at `n_chains`.

- seed:

  Random seed for reproducibility. Per-chain seeds are based on this
  seed, so a full multi-chain run is reproducible.

- verbose:

  Whether to print progress information (default: TRUE).

- log_performance:

  Whether to log performance metrics (default: FALSE).

- store_mcmc:

  Whether to store full MCMC samples (default: FALSE). Only post-burn-in
  iterations are stored.

- ...:

  Additional model-specific parameters.

- e1:

  Error parameter for categorical model (default: 0.05)

- e2:

  Error parameter for categorical model (default: 0.05)

- dict_prior:

  Dictionary prior for the multinomial model; see `model`.

## Value

An object of class `snp_slice_results` containing:

- `chains`: Per-chain results, all stored the same way. Each holds that
  chain's MAP estimate (`map_allocation_matrix` (A),
  `map_dictionary_matrix` (D)) and final-sample estimate
  (`final_allocation_matrix`, `final_dictionary_matrix`), plus
  `mcmc_samples` (if store_mcmc = TRUE), `diagnostics`, `convergence`,
  and `seed`

- `best_chain`: Index of the chain with the highest MAP log posterior

- `parameters`: MCMC settings used

- `model_info`: Model specification

The object holds no estimates of its own. Reach a chain's estimates with
[`get_chain()`](https://plasmogenepi.github.io/snp.slicer/reference/get_chain.md),
which defaults to the best chain, or with
[`extract_allocations()`](https://plasmogenepi.github.io/snp.slicer/reference/extract_allocations.md)
/
[`extract_strains()`](https://plasmogenepi.github.io/snp.slicer/reference/extract_strains.md);
every diagnostic function also takes a `chain` argument.
[`compare_chains()`](https://plasmogenepi.github.io/snp.slicer/reference/compare_chains.md)
summarises all chains.

## The multinomial model

The other models describe a target with a single number per specimen:
the dictionary is binary, so `(A %*% D) / rowSums(A)` gives the fraction
of a specimen's strains carrying the alternate allele, and the reference
allele is whatever is left over. That only works for two alleles, so
targets with more are dropped.

The multinomial model instead describes a target with a probability
*vector*, one entry per allele. The dictionary holds an allele code per
strain per target rather than a bit; expanding it to one indicator
column per (target, allele) pair makes each entry of
`(A %*% Dx) / rowSums(A)` the fraction of a specimen's strains carrying
that allele, and the entries for one target sum to one. The observed
read counts across a target's alleles are then multinomial in that
vector. Nothing limits the vector's length, so every allele at every
target is kept.

With exactly two alleles the vector is `(1 - p, p)` and the kernel
reduces to the binomial one, so the two models agree on biallelic data.

Parameters read from `...`: `dict_prior` sets the per-target categorical
prior on dictionary entries, either `"empirical"` (default; pooled
allele read fractions with a pseudocount of one), `"uniform"`, or a list
of per-target probability vectors. `rho`, used by the count models, does
not apply.

Implementation: the allocation update runs in compiled code as for the
other models, while the dictionary update is an R Gibbs step over
alleles.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with read count data
data <- list(
  read1 = matrix(c(10, 5, 15, 8), nrow = 2),
  read0 = matrix(c(90, 95, 85, 92), nrow = 2)
)

result <- snp_slice(data, model = "negative_binomial", n_sample = 1000)

# Extract results
strains <- extract_strains(result)
allocations <- extract_allocations(result)
} # }
```
