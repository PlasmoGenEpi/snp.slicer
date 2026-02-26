# Update single element of allocation matrix A

Update single element of allocation matrix A

## Usage

``` r
slice_update_a_local(
  state,
  i,
  k,
  model_obj,
  full_row_i = NULL,
  row_sum_i = NULL,
  kstar = NULL
)
```

## Arguments

- state:

  Current state

- i:

  Host index

- k:

  Strain index

- model_obj:

  Model object

- full_row_i:

  Optional precomputed row contribution `A[i,]` %\*% D (from
  slice_update_a)

- row_sum_i:

  Optional precomputed sum of `A[i,]` (from slice_update_a)

- kstar:

  Optional current last-active-feature index (from slice_update_a)

## Value

Updated state
