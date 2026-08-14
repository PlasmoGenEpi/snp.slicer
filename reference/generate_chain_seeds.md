# Draw independent per-chain seeds

Generate a separate random seed for each chain, based on the caller's
random number generator to maintain reproducibility.

## Usage

``` r
generate_chain_seeds(n_chains)
```

## Arguments

- n_chains:

  Number of chains

## Value

Integer vector of length `n_chains`

## Details

Seeds are drawn from the caller's RNG stream, so setting `seed` in
[`snp_slice`](https://plasmogenepi.github.io/snp.slicer/reference/snp_slice.md)
makes the whole multi-chain run reproducible while each chain still
starts from a different point. Because each chain receives an explicit
integer seed rather than inheriting global RNG state, the same seeding
scheme works for a compiled kernel with its own RNG.
