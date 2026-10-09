# R reference kernel adapter

R reference kernel adapter

## Usage

``` r
mcmc_kernel_r(model_obj = NULL)
```

## Arguments

- model_obj:

  Optional model object. The multinomial model has its own allocation
  and dictionary updates (integer dictionary, Gibbs over alleles); every
  other model, and a `NULL`, gets the binary reference updates.
