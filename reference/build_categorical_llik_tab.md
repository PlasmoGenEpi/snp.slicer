# Categorical Model for SNP-Slice

Implementation of the categorical observation model for SNP-Slice. This
model handles categorical observations (0, 0.5, 1) with error
parameters. Build categorical log-likelihood lookup table

## Usage

``` r
build_categorical_llik_tab(e1, e2)
```

## Arguments

- e1:

  False-positive error rate

- e2:

  False-negative / mixture error rate

## Value

3x3 matrix: rows = proportion bins, cols = observation bins
