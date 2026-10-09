# Expand an integer dictionary to one-hot allele-slot columns

Expand an integer dictionary to one-hot allele-slot columns

## Usage

``` r
expand_dictionary(D, model_obj)
```

## Arguments

- D:

  Integer dictionary, strains in rows, one allele code per target.

- model_obj:

  Multinomial model object carrying the expanded layout.

## Value

Binary matrix with `nrow(D)` rows and `model_obj$L` columns; row `k` has
a 1 in the slot of the allele strain `k` carries at each target.
