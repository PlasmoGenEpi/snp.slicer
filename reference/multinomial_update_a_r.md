# Allocation update for the multinomial model (R kernel)

Reuses
[`slice_update_a_r`](https://plasmogenepi.github.io/snp.slicer/reference/slice_update_a_r.md)
on the expanded layout: the integer dictionary is swapped for its
one-hot expansion and the observation matrix for the expanded allele
counts, so the generic per-cell Metropolis step evaluates the
multinomial likelihood unchanged.

## Usage

``` r
multinomial_update_a_r(state, model_obj)
```

## Arguments

- state:

  Current state

- model_obj:

  Model object

## Value

Updated state
