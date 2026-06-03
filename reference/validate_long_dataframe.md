# Structural validation for long-format data.frames

Checks shared by every long-format input (categorical and count models):
each allele must appear at most once per specimen and target, and at
least one target must be biallelic (exactly two observed alleles) so the
model has variation to resolve strains. Monomorphic targets are allowed
alongside biallelic ones; targets with more than two alleles are dropped
downstream by
[`load_dataframe`](https://plasmogenepi.github.io/snp.slicer/reference/load_dataframe.md).
The caller
([`validate_input_data`](https://plasmogenepi.github.io/snp.slicer/reference/validate_input_data.md))
guarantees the key columns are present before this is called.

## Usage

``` r
validate_long_dataframe(data, params)
```

## Arguments

- data:

  Long-format data.frame.

- params:

  List of resolved column names with elements `specimen_id_col`,
  `target_id_col`, `target_value_col`.

## Value

Invisibly `TRUE`; otherwise throws an error.
