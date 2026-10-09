# Log prior of the dictionary under the model's prior

Dispatches to `model_obj$logprior_d` when the model defines one and to
the Bernoulli(`rho`) prior in
[`logprior_d`](https://plasmogenepi.github.io/snp.slicer/reference/logprior_d.md)
otherwise.

## Usage

``` r
dictionary_logprior(D, model_obj)
```

## Arguments

- D:

  Dictionary matrix

- model_obj:

  Model object

## Value

Log prior probability
