# Resolve exceptions for multinomial model

An exception is a specimen with reads for an allele that none of its
assigned strains carries, which gives the likelihood `-Inf`. Each is
repaired by assigning the specimen an existing strain carrying that
allele at that target (preferring updatable strains, `k >= kmin`) or, if
there is none, appending a new strain equal to the specimen's most-read
genotype with that one allele substituted.

## Usage

``` r
multinomial_resolve_exceptions(state, model_obj)
```

## Arguments

- state:

  Current state

- model_obj:

  Model object

## Value

Updated state
