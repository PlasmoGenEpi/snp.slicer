# Compare chains from a multi-chain run

Compare chains from a multi-chain run

## Usage

``` r
compare_chains(results)
```

## Arguments

- results:

  SNP-Slice results object

## Value

A data frame with one row per chain: chain, seed, iterations_run,
map_logpost, final_logpost, map_kstar, n_strains, gap_converged, and
whether the chain was selected as best

## Examples

``` r
result <- load_example_results()
compare_chains(result)
#>   chain       seed iterations_run map_logpost final_logpost map_kstar n_strains
#> 1     1 1235143119           1250   -74834.45     -74968.33       113        52
#> 2     2 1756553742           1250   -75096.07     -75202.88       113        48
#> 3     3 1891765162           1250   -74918.26     -75013.13       113        53
#>   gap_converged  best
#> 1         FALSE  TRUE
#> 2         FALSE FALSE
#> 3         FALSE FALSE
```
