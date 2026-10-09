# Multinomial Model for SNP-Slice

Observation model for targets with any number of alleles. Each strain
carries exactly one allele per target, stored in the dictionary as an
integer allele code (`0` for slot 1, `1` for slot 2, and so on). For a
specimen carrying a set of strains, the expected allele proportions at a
target are the fractions of those strains carrying each allele, and the
observed read counts across alleles are multinomial with those
proportions.

With two alleles per target this is exactly the binomial model, so the
multinomial model nests it. The dictionary prior is a per-target
categorical distribution over alleles rather than the Bernoulli(`rho`)
prior of the biallelic models.

The allocation update, which dominates run time, runs in compiled code
via
[`mcmc_kernel_multinomial_cpp`](https://plasmogenepi.github.io/snp.slicer/reference/mcmc_kernel_multinomial_cpp.md)
when the compiled kernel is available;
[`multinomial_update_a_r`](https://plasmogenepi.github.io/snp.slicer/reference/multinomial_update_a_r.md)
is the R reference. The dictionary update is always the R Gibbs step
[`multinomial_update_d_r`](https://plasmogenepi.github.io/snp.slicer/reference/multinomial_update_d_r.md).

## Expanded layout

Internally, per-target allele counts are stored side by side in one
`N x L` matrix (`counts_exp`), where `L` is the total number of allele
slots across targets. `col_offset[p]` is the number of columns before
target `p`'s block, `col_locus` maps each column to its target, and
`col_allele` gives the 0-based allele code of each column.
[`expand_dictionary()`](https://plasmogenepi.github.io/snp.slicer/reference/expand_dictionary.md)
turns an integer dictionary into the matching one-hot `K x L` matrix, so
`A %*% Dx` counts, per specimen and allele slot, the strains carrying
that allele.
