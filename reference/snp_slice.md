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
  n_mcmc = 10000,
  burnin = NULL,
  alpha = 2.6,
  rho = 0.5,
  threshold = 0.001,
  gap = NULL,
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
```

## Arguments

- data:

  Input data. Can be a matrix, data.frame, or file path. For read count
  data, should be a list with elements `read1` and `read0` (or `total`).
  For categorical data, should be a matrix with values 0, 0.5, or 1.

- model:

  Observation model to use. Options: "categorical", "poisson",
  "binomial", "negative_binomial" (default).

- n_mcmc:

  Number of MCMC iterations (default: 10000).

- burnin:

  Burn-in period. If NULL, defaults to n_mcmc/2.

- alpha:

  IBP concentration parameter (default: 2.6).

- rho:

  Dictionary sparsity parameter (default: 0.5).

- threshold:

  Threshold for identifying single infections (default: 0.001).

- gap:

  Early stopping threshold. If NULL, runs for full n_mcmc iterations.

- seed:

  Random seed for reproducibility.

- verbose:

  Whether to print progress information (default: TRUE).

- log_performance:

  Whether to log performance metrics (default: FALSE).

- store_mcmc:

  Whether to store full MCMC samples (default: FALSE).

- ...:

  Additional model-specific parameters.

- e1:

  Error parameter for categorical model (default: 0.05)

- e2:

  Error parameter for categorical model (default: 0.05)

## Value

An object of class `snp_slice_results` containing:

- `allocation_matrix`: Binary allocation matrix (A)

- `dictionary_matrix`: Binary dictionary matrix (D)

- `mcmc_samples`: MCMC samples (if store_mcmc = TRUE)

- `diagnostics`: Convergence diagnostics

- `parameters`: Model parameters used

- `model_info`: Model specification

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with read count data
data <- list(
  read1 = matrix(c(10, 5, 15, 8), nrow = 2),
  read0 = matrix(c(90, 95, 85, 92), nrow = 2)
)

result <- snp_slice(data, model = "negative_binomial", n_mcmc = 1000)

# Extract results
strains <- extract_strains(result)
allocations <- extract_allocations(result)
} # }
```
