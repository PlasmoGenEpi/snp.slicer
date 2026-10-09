# Vectorised form of [`categorical_prop_bin()`](https://plasmogenepi.github.io/snp.slicer/reference/categorical_prop_bin.md)

The exception detector in
[`categorical_resolve_exceptions()`](https://plasmogenepi.github.io/snp.slicer/reference/categorical_resolve_exceptions.md)
must classify proportions exactly as the likelihood does, so both call
into the same bin edges. Keeping a separate scalar and matrix form is
deliberate: the scalar one is called per-cell in the likelihood's inner
loop.

## Usage

``` r
categorical_prop_bin_vec(prop)
```

## Arguments

- prop:

  Numeric vector or matrix of proportions.

## Value

Integer object of the same shape, with values 1, 2 or 3.
