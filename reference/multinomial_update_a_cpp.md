# Compiled allocation update for the multinomial model

Calls the locus-aware compiled update: for each specimen and strain, the
likelihood ratio between carrying and not carrying the strain involves
only the allele slot the strain carries at each target, so the compiled
loop visits one slot per target and takes its log terms from a table. It
reads the integer dictionary directly; no expanded layout is built.

## Usage

``` r
multinomial_update_a_cpp(state, model_obj)
```

## Arguments

- state:

  Current state

- model_obj:

  Model object

## Value

Updated state
