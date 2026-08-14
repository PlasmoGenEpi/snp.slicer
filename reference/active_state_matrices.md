# Extract allocation and dictionary matrices from a chain state

Extract allocation and dictionary matrices from a chain state

## Usage

``` r
active_state_matrices(state)
```

## Arguments

- state:

  A chain state (e.g. the MAP state or the final sample), carrying full
  `A` and `D` matrices

## Value

List with `allocation_matrix` and `dictionary_matrix`, restricted to
strains with at least one host
