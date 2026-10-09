# Resolve kernel adapter and attach cpp obs cache on `model_obj`

Call before a manual
[`slice_iter()`](https://plasmogenepi.github.io/snp.slicer/reference/slice_iter.md)
loop when not using
[`run_chain()`](https://plasmogenepi.github.io/snp.slicer/reference/run_chain.md).
[`run_chain()`](https://plasmogenepi.github.io/snp.slicer/reference/run_chain.md)
and
[`snp_slice()`](https://plasmogenepi.github.io/snp.slicer/reference/snp_slice.md)
invoke this automatically.

## Usage

``` r
setup_mcmc_kernel(
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

Updated model object (same list, with `kernel` and cache fields set).
