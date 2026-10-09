# Calculate allele counts for a single sample

Calculate allele counts for a single sample

## Usage

``` r
calculate_allele_counts_single(
  A,
  D,
  snp_indices,
  r0_values,
  r1_values,
  sep = "|",
  model = NULL,
  allele_labels = NULL,
  y_slot = 1L
)
```

## Arguments

- allele_labels:

  Optional list, one character vector per entry of `snp_indices`, giving
  the label of each allele code at that target (multinomial model). When
  supplied it overrides `r0_values` and `r1_values`.

- y_slot:

  For count models, the allele slot whose counts were modelled as `y`: 1
  (long-format input, `y = read0`) or 2 (read0/read1 list input,
  `y = read1`). A dictionary entry of 1 always means the strain carries
  the `y` allele, so this decides which label that is.
