# Load a long-format data.frame keeping every allele at every target

Multi-allelic counterpart of
[`load_dataframe`](https://plasmogenepi.github.io/snp.slicer/reference/load_dataframe.md)
for the multinomial model. No target is dropped for having more than two
alleles. Within each target, alleles are ordered by descending total
`target_count` across specimens (ties broken by `target_value`), so slot
1 is the population major allele and the dictionary code `0` denotes it.

## Usage

``` r
load_dataframe_multiallelic(
  data,
  model = "multinomial",
  target_id_col = "target_id",
  target_value_col = "target_value",
  specimen_id_col = "specimen_id",
  target_count_col = "target_count",
  ...
)
```

## Arguments

- data:

  Input dataframe with columns:

  specimen_id

  :   Specimen ID

  target_id

  :   Target ID

  target_value

  :   Target value (allele)

  target_count

  :   Target count

  Input is assumed to have passed
  [`validate_input_data`](https://plasmogenepi.github.io/snp.slicer/reference/validate_input_data.md):
  the four columns are present, no `(specimen, target, value)` row is
  duplicated, and at least one target is biallelic.

- model:

  Model type

- target_id_col:

  Name of the target ID column

- target_value_col:

  Name of the target value column

- specimen_id_col:

  Name of the specimen ID column

- target_count_col:

  Name of the target count column

- ...:

  Ignored; absorbs model parameters such as `dict_prior` that travel
  alongside column-name overrides.

## Value

Processed data list; see
[`build_multiallelic_processed`](https://plasmogenepi.github.io/snp.slicer/reference/build_multiallelic_processed.md)
for the fields.
