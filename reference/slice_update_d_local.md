# Update single element of dictionary matrix D

Update single element of dictionary matrix D

## Usage

``` r
slice_update_d_local(state, k, p, model_obj, an = NULL, ad0 = NULL)
```

## Arguments

- state:

  Current state

- k:

  Strain index

- p:

  SNP index

- model_obj:

  Model object

- an:

  Optional precomputed row sums of A (from slice_update_d)

- ad0:

  Optional precomputed leave-k-out contribution for this (k,p)

## Value

Updated state
