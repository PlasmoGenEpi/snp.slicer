# Load data from a dataframe

Load data from a dataframe

## Usage

``` r
load_dataframe(
  data,
  model,
  unknown_target_value = "?",
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

  :   Target value

  target_count

  :   Target count (optional)

  For each target, there are at exactly 2 taget values observed. If
  there is only one, the second value is set to unknown_target_value.
  Target count is required if model is not "categorical".

- model:

  Model type

- unknown_target_value:

  Value to use for unknown targets

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
