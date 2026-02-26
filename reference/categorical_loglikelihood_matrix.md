# Categorical Model for SNP-Slice

Implementation of the categorical observation model for SNP-Slice. This
model handles categorical observations (0, 0.5, 1) with error
parameters. Log-likelihood for categorical model (matrix version)

## Usage

``` r
categorical_loglikelihood_matrix(A, D, model_obj)
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
