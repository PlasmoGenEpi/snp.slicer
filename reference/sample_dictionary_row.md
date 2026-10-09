# Draw one dictionary row from the model's prior

Models that define `sample_dict_row` (the multinomial model) draw allele
codes from their per-target priors; every other model draws independent
Bernoulli(`rho`) bits.

## Usage

``` r
sample_dictionary_row(model_obj)
```

## Arguments

- model_obj:

  Model object

## Value

Numeric vector of length `model_obj$P`
