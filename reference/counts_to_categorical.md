# Convert read count matrices to categorical observation matrix

Maps reference (read0) and alternate (read1) counts to categorical
states: 0 (ref only), 1 (alt only), 0.5 (both), or NA (zero total).

## Usage

``` r
counts_to_categorical(read0, read1)
```

## Arguments

- read0:

  Matrix of reference allele counts (specimens x targets).

- read1:

  Matrix of alternate allele counts (same dimensions as read0).

## Value

Numeric matrix of same shape with values in {0, 0.5, 1, NA}, dimnames
preserved.
