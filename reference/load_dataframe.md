# Load data from a dataframe

Load data from a dataframe

## Usage

``` r
load_dataframe(
  data,
  model,
  target_id_col = "target_id",
  target_value_col = "target_value",
  specimen_id_col = "specimen_id",
  target_count_col = "target_count"
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

## Value

Processed data list with y, r, and metadata

## Details

The alleles for each target are sorted and given indices such that
target_idx 1 (read0) is the minor allele and target_idx 2 (read1) the
major. For monomorphic targets, with only one observed allele, both
slots carry that same label and read1 is zero (because no second allele
was observed). A specimen with no reads at a target (i.e., total_count
0) is encoded as `NA` for both slots, thus distinguishing a missing
genotype from a homozygous call.

When `model == "categorical"`, long-format data.frames are handled by
`load_dataframe_categorical`, which builds read0/read1 matrices from the
same layout and then converts counts to categorical observations. For
categorical data the allele order is determined by descending
`target_value` (not by count), so the lexicographically larger allele
label occupies target_idx 1 (read0). A specimen observed with only that
allele yields `y = 0`, with only the other allele `y = 1`, with both
`y = 0.5`, and with no reads `y = NA`.
