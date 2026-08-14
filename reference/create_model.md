# Create model object

Create model object

## Usage

``` r
create_model(
  model,
  processed_data,
  alpha = 2.6,
  rho = if (model == "categorical") 0.5 else NULL,
  ...
)
```

## Arguments

- model:

  Model type

- processed_data:

  Processed data

- ...:

  Additional model parameters

## Value

Model object with initialization and likelihood functions
