# De-duplicate rows of a matrix of allele codes

Like
[`remove_duplicates`](https://plasmogenepi.github.io/snp.slicer/reference/remove_duplicates.md)
but keyed on exact row equality, so it works for integer allele codes as
well as binary patterns.

## Usage

``` r
remove_duplicate_rows(dd)
```

## Arguments

- dd:

  Numeric matrix, one candidate strain per row

## Value

List with `assignments` (row -\> unique-row index) and `D` (the unique
rows, in order of first appearance)
