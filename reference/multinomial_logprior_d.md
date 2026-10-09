# Log prior of an integer dictionary under per-target categorical priors

Log prior of an integer dictionary under per-target categorical priors

## Usage

``` r
multinomial_logprior_d(D, model_obj)
```

## Arguments

- D:

  Integer dictionary matrix

- model_obj:

  Model object with `log_dict_prior_pad`, a `max(n_alleles) x P` matrix
  of log prior probabilities.

## Value

Log prior probability
