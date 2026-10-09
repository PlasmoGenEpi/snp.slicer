# Strains a host could plausibly be carrying, given its own calls

A mixed host used to be initialized carrying every strain in the
dictionary. That is not a biological state – on a 400-specimen panel it
starts hosts at a COI of ~289 – and it manufactures likelihood
exceptions that cannot be repaired. `llik_tab` is `-Inf` at (prop bin 3,
y == 0), bin 3 being `prop > 0.99`; a host carrying all K strains at a
locus where K-1 of them carry the alternate sits at `(K-1)/K`, which is
inside bin 3 once K \> 100.
[`categorical_resolve_exceptions()`](https://plasmogenepi.github.io/snp.slicer/reference/categorical_resolve_exceptions.md)
repairs such a cell by adding a single reference-carrying strain, which
only reaches bin 2 when K \<= 100, so above that the repair is
arithmetically incapable of clearing the exception and initialization
aborts at the iteration cap.

## Usage

``` r
categorical_compatible_strains(y_row, D, fallback)
```

## Arguments

- y_row:

  Observed calls for one host (0, 1, 0.5 or NA per locus).

- D:

  Strain dictionary, strains in rows.

- fallback:

  Strain index to use if nothing matches.

## Value

Integer vector of strain indices.

## Details

Restricting a host to the strains that match its own homozygous calls
fixes this at the source: at every locus called 0 or 1 the host then
carries only strains with that allele, so the proportion is exactly 0 or
1 and lands in the bin the observation agrees with. Heterozygous and
missing loci are left unconstrained, since any mixture is consistent
with them.

The result is never empty: the dictionary is built by de-duplicating the
hosts' own resolved genotypes, so a host's own strain always matches it
at every homozygous locus. `fallback` guards the degenerate case anyway.
