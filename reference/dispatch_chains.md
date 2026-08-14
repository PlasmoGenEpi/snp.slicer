# Dispatch a per-chain function, possibly in parallel

Run the provided function on each chain, in parallel if `n_cores` \> 1.

## Usage

``` r
dispatch_chains(n_chains, n_cores, fun)
```

## Arguments

- n_chains:

  Number of chains

- n_cores:

  Number of cores to use

- fun:

  Function of a single argument (the chain index)

## Value

List of chain results

## Details

Uses forking
([`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html)) where
available and a PSOCK cluster on Windows. Falls back to `lapply` when a
single core is requested.
