# Assemble processed multi-allelic data from per-target count matrices

Assemble processed multi-allelic data from per-target count matrices

## Usage

``` r
build_multiallelic_processed(
  counts_list,
  allele_labels,
  specimen_ids,
  target_ids,
  model
)
```

## Arguments

- counts_list:

  List with one `N x M_p` count matrix per target, columns in
  allele-slot order.

- allele_labels:

  List of character vectors, the allele label of each slot per target.

- specimen_ids, target_ids:

  Row and target identifiers.

- model:

  Model name recorded in the output.

## Value

A processed data list with the fields the biallelic loaders produce (`y`
holds slot-1 counts, `r` the per-target totals, `NA` where a specimen
has no reads at a target) plus the expanded layout: `counts_exp`
(`N x L`), `n_alleles`, `col_offset`, `col_locus`, `col_allele`, and
`allele_labels`. `r0_values`/`r1_values` carry the slot-1 and slot-2
labels for code that expects them; slot 2 repeats slot 1 at monomorphic
targets.
