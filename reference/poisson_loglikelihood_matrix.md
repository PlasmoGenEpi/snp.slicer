# Poisson Model for SNP-Slice

Implementation of the Poisson observation model for SNP-Slice. This
model assumes Y ~ Poisson(R × proportions). Log-likelihood for Poisson
model (matrix version)

## Usage

``` r
poisson_loglikelihood_matrix(A, D, model_obj)
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
