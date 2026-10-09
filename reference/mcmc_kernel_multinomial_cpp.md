# Compiled kernel adapter for the multinomial model

Like the biallelic compiled adapter, every update runs in compiled code
and `update_iter` fuses the slice-variable, allocation, dictionary and
stick-breaking steps into one call. The adapter also supplies `loglik`,
the compiled kernel term plus the model's multinomial coefficients.

## Usage

``` r
mcmc_kernel_multinomial_cpp(model_obj, fused = TRUE)
```

## Arguments

- model_obj:

  Multinomial model object

- fused:

  Use the fused compiled iteration (default). With `FALSE` the adapter
  runs the compiled allocation and dictionary updates with the R
  slice-variable and stick-breaking updates, which is the configuration
  comparable step for step with the R reference.

## Value

Kernel adapter list

## Details

The adapter's name is `"cpp_multinomial"` rather than `"cpp"` which is
the biallelic adaptor name:
[`slice_iter()`](https://plasmogenepi.github.io/snp.slicer/reference/slice_iter.md)
must not route it through the generic observation cache or the generic
total likelihood but instead via the multinomial one.
