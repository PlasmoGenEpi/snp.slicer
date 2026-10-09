# Build the multinomial model object

Build the multinomial model object

## Usage

``` r
build_multinomial_model(processed_data, alpha, dict_prior = "empirical")
```

## Arguments

- processed_data:

  Output of
  [`load_dataframe_multiallelic`](https://plasmogenepi.github.io/snp.slicer/reference/load_dataframe_multiallelic.md)
  or the list-input path of
  [`preprocess_data`](https://plasmogenepi.github.io/snp.slicer/reference/preprocess_data.md)
  with `model = "multinomial"`.

- alpha:

  IBP concentration parameter

- dict_prior:

  See
  [`resolve_dict_prior`](https://plasmogenepi.github.io/snp.slicer/reference/resolve_dict_prior.md).

## Value

Model object list (class assigned by
[`create_model`](https://plasmogenepi.github.io/snp.slicer/reference/create_model.md)).
