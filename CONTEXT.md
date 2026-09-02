# SNP-Slice (snp.slicer)

Bayesian nonparametric resolution of multi-strain infections from SNP/target data. This glossary covers domain terms for data representation and input — not implementation details.

## Language

**Target**:
A genotyped site (SNP or equivalent) identified by `target_id`. Used in long-format input.
_Avoid_: Locus (unless speaking informally), SNP (when the input is not strictly biallelic SNPs)

**Specimen**:
A biological sample from which reads were obtained, identified by `specimen_id`.
_Avoid_: Host, individual (in user-facing docs; "individual" may appear in internal MCMC state)

**Allele slot 1 / allele slot 2**:
The two fixed columns in internal read-count matrices. Slot 1 counts are stored in `y` (read0); slot 2 counts are derived as `r - y` (read1). Allele labels for each slot are `r0_values` and `r1_values`.
_Avoid_: Overloading "reference/alternate" unless slot assignment is explicitly reference-first

**Allele slot ordering (long-format input)**:
Rules for which allele occupies slot 1 vs slot 2 when the package assigns order:
- Count models (Poisson, binomial, negative binomial): slot 1 = population **minor** allele (lowest total `target_count` across specimens); slot 2 = major allele; lexicographic order breaks ties. Order is global per target (not per specimen).
- Categorical model: slot 1 = lexicographically **larger** allele; slot 2 = the other.
- Matrix input (`list(read0, read1)`): no reordering; user-supplied orientation is preserved.

**Reference / alternate (categorical encoding only)**:
In the categorical observation model, slot 1 allele = reference (observation `0`), slot 2 allele = alternate (observation `1`), both present = `0.5`. These names apply to categorical encoding semantics, not to slot assignment in count models.
_Avoid_: Using "ref/alt" to describe count-model minor/major slot assignment

**Monomorphic target**:
A target at which only one allele is observed across the entire dataset. It constrains every strain to carry that allele (homozygosity at an invariant site) but provides no allelic variation for strain discrimination.
_Avoid_: Calling it "uninformative" (it is informative as a constraint, just not for separating strains)

**Monomorphic target encoding**:
Both allele slots carry the same label; all reads go in slot 1 and slot 2 is a structural zero (count 0, not missing). Specimens with zero total reads at the target are encoded as `NA` (missing genotype), not as homozygous.
_Avoid_: Using a placeholder allele (e.g. `"?"`) for the second slot

**Monomorphic target (categorical model)**:
An observed monomorphic call with reads encodes as `0` (slot 1 present, slot 2 absent). Because both slots share the same allele label, this is equivalent to homozygous for the sole allele — not biallelic ref-only in the usual sense.
_Avoid_: Encoding as `0.5` or dropping monomorphic targets from categorical input

**Biallelic target**:
A target with exactly two distinct alleles observed in the dataset. At least one biallelic target is required for inference; a dataset containing only monomorphic targets is rejected.
_Avoid_: Requiring every target to be biallelic

**Rho (`rho`)**:
Dictionary sparsity prior — per-SNP probability that a strain carries the alternate allele in the dictionary matrix `D`.
_Avoid_: Confusing with correlation or rate parameters

**Default rho**:
- Categorical model: fixed `0.5`.
- Count models (Poisson, binomial, negative binomial): `"MAF"` — resolved to the global minor-allele read fraction **at biallelic targets only** (`sum(y) / sum(r)` excluding monomorphic targets). Monomorphic coverage must not enter the ratio. At least one biallelic target must have positive total reads (validated at input). If MAF is still undefined at model creation, fall back to `0.5` with a warning.
_Avoid_: Including monomorphic targets in the MAF ratio (they inflate `rho` toward 1)

**Polyallelic target**:
A target with more than two distinct alleles observed in the dataset. SNP-Slice is biallelic-only; such targets are dropped during loading with a warning listing the dropped `target_id`s.
_Avoid_: Silently discarding them, or erroring on the entire dataset for a single triallelic site

**Long-format row uniqueness**:
At most one row per `(specimen_id, target_id, target_value)`. Duplicate keys are invalid input and must error — never silently aggregated.
_Avoid_: Summing duplicate rows during load

**Missing genotype**:
A `(specimen, target)` pair with no reads after grid completion (zero total coverage across both allele slots). Encoded as `NA` in internal matrices. Sparse long-format input — rows only where reads exist — is valid; absent pairs are completed and treated as missing.
_Avoid_: Treating absence as a structural zero observation, or requiring dense long-format input

## Example dialogue


**Dev:** For this long-format dataframe with NB model, which allele goes in slot 1?
**Expert:** Sum counts across all specimens per allele at that target. Whichever allele has fewer total reads is slot 1 — that's the minor allele SNP-Slice expects in the first column.
**Dev:** Same target but categorical model?
**Expert:** Ignore counts for ordering. Lexicographically larger allele is slot 1; the other is slot 2. Counts still determine whether the observation is 0, 0.5, 1, or NA.
**Dev:** User passes read0/read1 matrices directly?
**Expert:** We don't reorder. Whatever they put in read0 is slot 1; read1 is slot 2.
**Dev:** Specimen s3 has no rows at all for target t2?
**Expert:** Sparse format is fine — we complete the grid. Zero total coverage means missing genotype, NA. That's different from a monomorphic G call where reads are present.
