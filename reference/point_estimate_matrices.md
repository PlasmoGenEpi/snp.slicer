# Access the allocation and dictionary matrices for a point estimate

Selects the matrices backing a point estimate, either the MAP or the
final sample of the chain.

## Usage

``` r
point_estimate_matrices(results, estimate)
```

## Arguments

- results:

  A single-chain `snp_slice_results` object (as returned by
  [`get_chain`](https://plasmogenepi.github.io/snp.slicer/reference/get_chain.md)).

- estimate:

  Either `"final_sample"` or `"map"`.

## Value

List with `A` (allocation matrix) and `D` (dictionary matrix).
