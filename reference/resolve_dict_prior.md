# Resolve per-target dictionary priors

Resolve per-target dictionary priors

## Usage

``` r
resolve_dict_prior(dict_prior, processed_data)
```

## Arguments

- dict_prior:

  `"empirical"` (pooled allele read fractions with a pseudocount of one
  per allele), `"uniform"`, or a list of numeric vectors, one per
  target, each of length `n_alleles[p]`.

- processed_data:

  Processed multi-allelic data.

## Value

List of probability vectors, one per target.
