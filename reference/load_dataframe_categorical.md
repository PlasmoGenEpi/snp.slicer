# Load long-format data.frame for categorical model

Uses the same join/completion logic as
[`load_dataframe`](https://plasmogenepi.github.io/snp.slicer/reference/load_dataframe.md)
to build read0/read1 matrices, then converts counts to categorical
observations (0, 0.5, 1, NA) via
[`counts_to_categorical`](https://plasmogenepi.github.io/snp.slicer/reference/counts_to_categorical.md).

## Usage

``` r
load_dataframe_categorical(
  data,
  model,
  target_id_col = "target_id",
  target_value_col = "target_value",
  specimen_id_col = "specimen_id",
  target_count_col = "target_count",
  ...
)
```

## Arguments

- data:

  Input dataframe (same column semantics as `load_dataframe`).

- model:

  Model type (must be `"categorical"`).

- target_id_col, target_value_col, specimen_id_col, target_count_col:

  Same as
  [`load_dataframe`](https://plasmogenepi.github.io/snp.slicer/reference/load_dataframe.md).

## Value

Processed data list with `data_type = "categorical"`, `y` the
categorical matrix, and `r = NULL`.
