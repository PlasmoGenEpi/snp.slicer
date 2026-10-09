# Log-likelihood for multinomial model (matrix version)

Includes the log multinomial coefficients, so values are comparable with
the binomial model's
[`dbinom()`](https://rdrr.io/r/stats/Binomial.html)-based likelihood on
two-allele data.

## Usage

``` r
multinomial_loglikelihood_matrix(A, D, model_obj)
```

## Arguments

- A:

  Allocation matrix

- D:

  Integer dictionary matrix

- model_obj:

  Model object containing data

## Value

Log-likelihood value
