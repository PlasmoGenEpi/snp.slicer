# Dictionary update for the multinomial model (R kernel)

Exact Gibbs step over each updatable dictionary cell: for strain `k` and
target `p`, the allele is redrawn from its full conditional given every
other cell, the current allocations, and the per-target prior. Only
specimens carrying strain `k` enter the likelihood ratio, and only the
slot each option adds carriers to (see
[`multinomial_cell_logweights`](https://plasmogenepi.github.io/snp.slicer/reference/multinomial_cell_logweights.md));
the specimen-by-slot carrier count matrix is kept incrementally. A cell
with a single feasible allele is set to it without a draw, and a strain
no specimen carries is redrawn from the prior.

## Usage

``` r
multinomial_update_d_r(state, model_obj)
```

## Arguments

- state:

  Current state

- model_obj:

  Model object

## Value

Updated state
