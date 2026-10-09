# Initialize state for multinomial model

Each specimen's most-read allele per target gives a candidate strain;
the de-duplicated candidates form the initial dictionary. A specimen is
a single infection when its top allele holds at least `1 - threshold` of
the reads at every observed target, and is assigned only its own strain.
A mixed specimen starts on every strain compatible with its own calls:
strains matching its allele at each target where it shows exactly one
allele. Alleles a specimen shows but no assigned strain carries are then
repaired by
[`multinomial_resolve_exceptions`](https://plasmogenepi.github.io/snp.slicer/reference/multinomial_resolve_exceptions.md).

## Usage

``` r
multinomial_initialize_state(model_obj, threshold = 0.001)
```

## Arguments

- model_obj:

  Model object

- threshold:

  Threshold for identifying single infections

## Value

Initialized state
