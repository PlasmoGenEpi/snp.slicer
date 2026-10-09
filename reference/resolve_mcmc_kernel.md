# Resolve MCMC kernel adapter for a model

Resolve MCMC kernel adapter for a model

## Usage

``` r
resolve_mcmc_kernel(
  model_obj,
  backend = getOption("snp.slicer.mcmc_kernel", "auto")
)
```

## Arguments

- model_obj:

  Model object from
  [`create_model()`](https://plasmogenepi.github.io/snp.slicer/reference/create_model.md).

- backend:

  One of `"auto"`, `"r"`, or `"cpp"`.

## Value

Kernel adapter list.
