# Log weights for the allele options of one dictionary cell

Pure function behind the Gibbs step in
[`multinomial_update_d_r`](https://plasmogenepi.github.io/snp.slicer/reference/multinomial_update_d_r.md).
Setting the cell to allele `m` adds the strain's carriers to slot `m`
only, so relative to a constant shared by every option the log
conditional of `m` is
`log_prior[m] + sum_h y[h, m] * (log(base[h, m] + 1) - log(base[h, m]))`.
The normalisation by each specimen's strain count is the same for every
option and drops out.

## Usage

``` r
multinomial_cell_logweights(base, yb, a_h, log_prior)
```

## Arguments

- base:

  Leave-strain-out carrier counts for the carrying specimens (specimens
  x alleles at this target).

- yb:

  Observed allele counts for those specimens at this target.

- a_h:

  The strain's allocation value for each of those specimens.

- log_prior:

  Log prior over the target's alleles.

## Value

Numeric vector of log weights, one per allele.

## Details

A slot with reads in some carrying specimen and no other carrier there
(`base == 0` with `y > 0`) makes every other option impossible: that
allele gets weight 0 on the log scale and the rest `-Inf`. In a valid
state only the current allele can be such a slot.
