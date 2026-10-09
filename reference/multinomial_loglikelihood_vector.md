# Log-likelihood for multinomial model (expanded-vector version)

Kernel term only: `sum(y * log(prop))` over observed cells with reads.
The multinomial coefficient is dropped because this function is used in
Metropolis and Gibbs ratios where it cancels. A cell with reads for an
allele no assigned strain carries contributes `-Inf`.

## Usage

``` r
multinomial_loglikelihood_vector(propvec, yvec, rvec = NULL)
```

## Arguments

- propvec:

  Expected allele proportions in the expanded layout (a vector for one
  specimen or a matrix for several).

- yvec:

  Observed allele counts in the same layout; `NA` marks a missing
  genotype.

- rvec:

  Ignored; kept for interface consistency with the other models.

## Value

Log-likelihood value
