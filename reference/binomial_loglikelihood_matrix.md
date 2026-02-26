# Binomial Model for SNP-Slice

Implementation of the binomial observation model for SNP-Slice. This
model assumes Y ~ Binomial(R, proportions). Log-likelihood for binomial
model (matrix version)

## Usage

``` r
binomial_loglikelihood_matrix(A, D, model_obj)
```

## Arguments

- A:

  Allocation matrix

- D:

  Dictionary matrix

- model_obj:

  Model object containing data

## Value

Log-likelihood value
