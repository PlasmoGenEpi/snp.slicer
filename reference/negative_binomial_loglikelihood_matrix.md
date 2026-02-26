# Negative Binomial Model for SNP-Slice

Implementation of the negative binomial observation model for SNP-Slice.
This is the recommended default model for most applications.
Log-likelihood for negative binomial model (matrix version)

## Usage

``` r
negative_binomial_loglikelihood_matrix(A, D, model_obj)
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
