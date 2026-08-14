# Introduction to SNP-Slice

## Overview

SNP-Slice is a Bayesian nonparametric method for resolving multi-strain
infections using slice sampling with stick-breaking construction. The
algorithm simultaneously unveils strain haplotypes and links them to
hosts from sequencing data.

This vignette demonstrates how to use SNP-Slice with the negative
binomial model using real example data.

## Installation

``` r

# Install from CRAN (when available)
install.packages("snp.slicer")

# Or install from GitHub
devtools::install_github("plasmogenepi/snp.slice")
```

## Quick start

Load the package, example data, and pre-computed results:

``` r

library(snp.slicer)

# Example data: read count matrices (hosts × SNPs)
data(example_snp_data, package = "snp.slicer")
read0 <- example_snp_data$read0
read1 <- example_snp_data$read1
cat("Data dimensions:", nrow(read0), "hosts ×", ncol(read0), "SNPs\n")
#> Data dimensions: 200 hosts × 96 SNPs

# Pre-computed results (negative binomial model, 3 chains)
result <- load_example_results()
cat("Pre-computed results loaded successfully!\n")
#> Pre-computed results loaded successfully!

# Basic inspection
print(result)
#> SNP-Slice Results
#> ================
#> Model: negative_binomial 
#> Dimensions: 200 hosts x 96 strains x 96 SNPs
#> Strains identified: 52 
#> Chains: 3 (best chain: 1 )
#> Gap Converged: No
summary(result)
#> SNP-Slice Results Summary
#> ========================
#> 
#> Model: negative_binomial 
#> Data dimensions: 200 hosts x 96 SNPs
#> Data type: read_counts 
#> 
#> Results:
#> - Number of strains identified: 52 
#> - Number of hosts: 200 
#> - Multiplicity of infection (MOI):
#>   - Mean MOI: 2.61 
#>   - Median MOI: 1 
#>   - Range: 1 - 11 
#>   - Single infections: 105 ( 52.5 %)
#>   - Mixed infections: 95 ( 47.5 %)
#> 
#> Chains: 3 (best chain: 1 )
#>  chain       seed iterations_run map_logpost final_logpost map_kstar n_strains
#>      1 1235143119           1250   -74834.45     -74968.33       113        52
#>      2 1756553742           1250   -75096.07     -75202.88       113        48
#>      3 1891765162           1250   -74918.26     -75013.13       113        53
#>  gap_converged  best
#>          FALSE  TRUE
#>          FALSE FALSE
#>          FALSE FALSE
#> 
#> Convergence:
#> - Iterations run: 1250 
#> - Samples retained (post-burn-in): 250 
#> - Gap Converged: No 
#> - Final log posterior: -74968.33 
#> - MAP log posterior: -74834.45 
#> - Final k*: 113 
#> - MAP k*: 113
```

## Understanding results

The main outputs are the **allocation matrix** (hosts × strains) and the
**dictionary matrix** (strains × SNPs). Use the extractors for
convenience:

``` r

# Strain and allocation summaries
strains <- extract_strains(result)
allocations <- extract_allocations(result)

cat("Strains identified:", strains$n_strains, "| SNPs:", strains$n_snps, "\n")
#> Strains identified: 52 | SNPs: 96
cat("Hosts:", allocations$n_hosts, "\n")
#> Hosts: 200
cat("Multiplicity of infection (MOI) summary:\n")
#> Multiplicity of infection (MOI) summary:
print(summary(allocations$multiplicity_of_infection))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    1.00    1.00    1.00    2.61    4.00   11.00

# Key matrices
A <- allocations$allocation_matrix   # Hosts × Strains
D <- strains$dictionary              # Strains × SNPs

cat("Allocation matrix:", dim(A), "| Dictionary matrix:", dim(D), "\n")
#> Allocation matrix: 200 52 | Dictionary matrix: 52 96
print("Sample allocation (first 3 hosts, first 5 strains):")
#> [1] "Sample allocation (first 3 hosts, first 5 strains):"
print(A[1:3, 1:5])
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    1    0    0    0    0
print("Sample dictionary (first 3 strains, first 8 SNPs):")
#> [1] "Sample dictionary (first 3 strains, first 8 SNPs):"
print(D[1:3, 1:8])
#>      site1 site2 site3 site4 site5 site6 site7 site8
#> [1,]     1     1     1     1     1     1     1     1
#> [2,]     1     1     1     1     1     1     0     1
#> [3,]     1     1     0     1     1     1     0     1
```

## Downstream analyses

### Posterior allele frequencies

Compute allele (or haplotype) frequencies for specific SNP sets. Use
`calculate_allele_frequencies` for a single set; use
`calculate_allele_frequencies_by_sets` for multiple sets at once. The
`estimate` argument chooses the estimator: `"final_sample"` (the final
sample of the chain, matching the SNP-Slice paper; the default) and
`"map"` (the maximum a posteriori allocation) are point estimates giving
literal counts (`allele`, `frequency`, `count`, `total_parasites`);
`"posterior"` gives a posterior mean frequency, SD, credible interval,
and per-sample mean count (`mean_count`, `n_samples`), so the table does
not depend on how many MCMC samples were used.

``` r

# Single target set: allele frequencies for SNPs 1, 5, and 10 (final sample)
allele_freqs <- calculate_allele_frequencies(result, c(1, 5, 10))
head(allele_freqs)
#>        allele  frequency count total_parasites
#> 8 ref|ref|ref 0.63428571   333             525
#> 2 ref|alt|alt 0.13523810    71             525
#> 3 alt|ref|alt 0.09142857    48             525
#> 4 ref|ref|alt 0.07428571    39             525
#> 7 alt|ref|ref 0.03619048    19             525
#> 6 ref|alt|ref 0.01714286     9             525

# Posterior: mean, SD, credible interval, and sample-size-invariant mean_count
if (!is.null(get_chain(result)$mcmc_samples)) {
  allele_freqs_post <- calculate_allele_frequencies(result, c(1, 5, 10), estimate = "posterior", n_samples = 50)
  head(allele_freqs_post)
}
#>                  allele  frequency frequency_sd frequency_lower frequency_upper
#> ref|ref|ref ref|ref|ref 0.63409753  0.002736118     0.631020945      0.63884162
#> ref|alt|alt ref|alt|alt 0.13803500  0.005609773     0.124605026      0.14621117
#> alt|ref|alt alt|ref|alt 0.08239172  0.009436347     0.066568441      0.09871358
#> ref|ref|alt ref|ref|alt 0.07315501  0.004944094     0.066568441      0.08405791
#> alt|ref|ref alt|ref|ref 0.04225196  0.007466272     0.032224515      0.04961832
#> ref|alt|ref ref|alt|ref 0.01446428  0.002814664     0.009527899      0.01915281
#>             mean_count n_samples
#> ref|ref|ref     332.32        50
#> ref|alt|alt      72.34        50
#> alt|ref|alt      43.18        50
#> ref|ref|alt      38.34        50
#> alt|ref|ref      22.14        50
#> ref|alt|ref       7.58        50

# Multiple target sets: pass a named list for one table per set
target_sets <- list(
  locus_a = c(1, 5),
  locus_b = c(10, 15),
  locus_c = c(20)
)
freqs_by_set <- calculate_allele_frequencies_by_sets(result, target_sets)
# Each element is a frequency table (point estimate: allele, frequency, count, total_parasites; posterior: adds frequency_sd, frequency_lower, frequency_upper, mean_count, n_samples)
print(freqs_by_set$locus_a)
#>    allele  frequency count total_parasites
#> 4 ref|ref 0.70857143   372             525
#> 2 ref|alt 0.15238095    80             525
#> 3 alt|ref 0.12761905    67             525
#> 1 alt|alt 0.01142857     6             525
```

### Individual COI with uncertainty

Per-host complexity of infection (COI) is the number of distinct strains
per host. Use
[`calculate_individual_coi()`](https://plasmogenepi.github.io/snp.slicer/reference/calculate_individual_coi.md)
for a point estimate (`estimate = "final_sample"` or `"map"`) or, when
MCMC samples were stored, `estimate = "posterior"` for the posterior
mean, SD, and a credible interval.

``` r

# Point estimate (from the final sample of the chain)
coi_final <- calculate_individual_coi(result, estimate = "final_sample")
head(coi_final)
#>   host_index    host_id coi_estimate coi_sd coi_lower coi_upper
#> 1          1 specimen_1            1     NA        NA        NA
#> 2          2 specimen_2            1     NA        NA        NA
#> 3          3 specimen_3            1     NA        NA        NA
#> 4          4 specimen_4            1     NA        NA        NA
#> 5          5 specimen_5            1     NA        NA        NA
#> 6          6 specimen_6            1     NA        NA        NA

# With posterior uncertainty (requires store_mcmc = TRUE)
if (!is.null(get_chain(result)$mcmc_samples)) {
  coi_post <- calculate_individual_coi(result, estimate = "posterior", n_samples = 50, interval = 0.95)
  head(coi_post)
  cat("\nSummary of COI posterior SD:\n")
  print(summary(coi_post$coi_sd))
}
#> 
#> Summary of COI posterior SD:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.00000 0.00000 0.04924 0.00000 0.99714
```

## Running Your Own Analysis

If you want to run the analysis yourself, here’s how to do it:

``` r

# Run SNP-Slice with the negative binomial model
data <- list(read0 = read0, read1 = read1)
result <- snp_slice(data, 
                   model = "negative_binomial",
                   n_sample = 2000,      # Number of retained (post-burn-in) samples
                   store_mcmc = TRUE,  # Store MCMC samples for diagnostics
                   verbose = TRUE)     # Show progress
```

## Parameter Tuning

You can adjust various parameters to optimize performance:

``` r

# Custom parameters
result_tuned <- snp_slice(data,
                         alpha = 1.5,        # IBP concentration parameter
                         rho = 0.3,          # Dictionary sparsity
                         threshold = 0.005,  # Single infection threshold
                         n_sample = 2000,
                         n_burnin = 500,       # Custom burn-in
                         gap = 10,           # Early stopping threshold
                         seed = 456,         # Reproducibility
                         verbose = FALSE)
```

## Convergence diagnostics

When you store MCMC samples (`store_mcmc = TRUE`),
[`convergence_diagnostics()`](https://plasmogenepi.github.io/snp.slicer/reference/convergence_diagnostics.md)
reports R-hat and effective sample size (ESS), pooled across all chains,
for permutation-invariant summaries of the run. R-hat is most
informative with several chains (run `snp_slice(..., n_chains = 3)`),
but a single chain still yields a split-R-hat.

``` r

if (!is.null(get_chain(result)$mcmc_samples)) {
  # R-hat and ESS for the whole run, one row per parameter
  print(convergence_diagnostics(result))

  # Trace plot of log posterior
  plot_convergence(result, type = "logpost")
} else {
  cat("MCMC samples not stored in results (set store_mcmc = TRUE when running snp_slice)\n")
}
#>    variable         mean    median         sd        q5       q95     rhat
#> 1   logpost -75055.92244 -75010.39 122.577982 -75254.45 -74902.59 2.187591
#> 2 n_strains     51.06133     52.00   2.157725     48.00     53.00 7.593054
#> 3     kstar    113.00000    113.00   0.000000    113.00    113.00       NA
#> 4    ktrunc    134.00000    134.00   0.000000    134.00    134.00       NA
#>   ess_bulk  ess_tail
#> 1 4.008605 41.678197
#> 2 3.322269  3.522944
#> 3       NA        NA
#> 4       NA        NA
```

![Trace plot of log posterior over MCMC
iterations.](introduction_files/figure-html/unnamed-chunk-9-1.png)

Constant parameters (e.g. `kstar`/`ktrunc` when the dictionary size has
settled) have zero variance, for which `posterior` returns `NaN` —
meaning “this parameter did not move”, not a failure.

Use `pars = "coi"` for per-host complexity of infection, which is
invariant to strain relabeling; element-wise diagnostics on the `A`/`D`
matrices are intentionally not provided because strain labels are not
identifiable.

[`as_draws_snp_slice()`](https://plasmogenepi.github.io/snp.slicer/reference/as_draws_snp_slice.md)
returns the underlying
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
if you want to apply other `posterior` functions directly. See
[`?convergence_diagnostics`](https://plasmogenepi.github.io/snp.slicer/reference/convergence_diagnostics.md)
and
[`?plot_convergence`](https://plasmogenepi.github.io/snp.slicer/reference/plot_convergence.md)
for options.

## Visualizing results

Let’s create some informative plots to better understand the analysis
results:

### Multiplicity of Infection (MOI) Distribution

``` r

# Get MOI from allocations (same as strain diversity)
allocations <- extract_allocations(result)
moi_df <- data.frame(
  moi = allocations$multiplicity_of_infection,
  host_id = 1:length(allocations$multiplicity_of_infection)
)

ggplot(moi_df, aes(x = moi)) +
  geom_histogram(binwidth = 1, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(aes(xintercept = mean(moi)), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of Multiplicity of Infection (MOI)",
    subtitle = "Number of Strains per Host",
    x = "Number of Strains per Host",
    y = "Number of Hosts",
    caption = paste("Mean MOI:", round(mean(moi_df$moi), 2), "| Max MOI:", max(moi_df$moi))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
```

![Histogram of multiplicity of infection (number of strains per host)
with mean indicated by a vertical
line.](introduction_files/figure-html/unnamed-chunk-10-1.png)

### Strain Frequency Analysis

``` r

# Calculate strain frequencies
strain_frequencies <- colSums(A)
freq_df <- data.frame(
  strain_id = 1:length(strain_frequencies),
  frequency = strain_frequencies
)

# Plot strain frequencies
ggplot(freq_df, aes(x = reorder(strain_id, frequency), y = frequency)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  labs(
    title = "Strain Frequencies Across Hosts",
    x = "Strain ID",
    y = "Number of Hosts Infected",
    caption = paste("Total strains:", length(strain_frequencies))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10)
  )
```

![Bar chart of strain frequencies across hosts (number of hosts infected
per strain).](introduction_files/figure-html/unnamed-chunk-11-1.png)

### Strain Pattern Heatmap

``` r

# Create a heatmap of the first 20 strains and first 30 SNPs
D_subset <- D[1:min(20, nrow(D)), 1:min(30, ncol(D))]

# Convert to long format for ggplot
heatmap_data <- expand.grid(
  strain = 1:nrow(D_subset),
  snp = 1:ncol(D_subset)
)
heatmap_data$value <- as.vector(D_subset)

ggplot(heatmap_data, aes(x = snp, y = strain, fill = factor(value))) +
  geom_tile() +
  scale_fill_manual(
    values = c("0" = "white", "1" = "darkblue"),
    labels = c("0" = "Reference", "1" = "Alternative"),
    name = "Allele"
  ) +
  labs(
    title = "Strain Patterns (First 20 Strains, First 30 SNPs)",
    x = "SNP Position",
    y = "Strain ID"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )
```

![Heatmap of strain genotypes: rows are strains, columns are SNPs; white
is reference allele, blue is
alternative.](introduction_files/figure-html/unnamed-chunk-12-1.png)

### Host-Strain Allocation Heatmap

``` r

# Create a heatmap of host-strain allocations (first 50 hosts, first 20 strains)
A_subset <- A[1:min(50, nrow(A)), 1:min(20, ncol(A))]

# Convert to long format for ggplot
allocation_data <- expand.grid(
  host = 1:nrow(A_subset),
  strain = 1:ncol(A_subset)
)
allocation_data$value <- as.vector(A_subset)

ggplot(allocation_data, aes(x = strain, y = host, fill = value)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white", 
    high = "red",
    name = "Allocation\nWeight"
  ) +
  labs(
    title = "Host-Strain Allocations (First 50 Hosts, First 20 Strains)",
    x = "Strain ID",
    y = "Host ID"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )
```

![Heatmap of host-strain allocation weights: rows are hosts, columns are
strains; color indicates allocation
strength.](introduction_files/figure-html/unnamed-chunk-13-1.png)

### Summary Statistics

``` r

# Create summary statistics
summary_stats <- data.frame(
  Metric = c(
    "Total Hosts",
    "Total SNPs", 
    "Total Strains",
    "Mean MOI",
    "Max MOI",
    "Single Infections",
    "Multiple Infections"
  ),
  Value = c(
    nrow(A),
    ncol(D),
    ncol(A),
    round(mean(allocations$multiplicity_of_infection), 2),
    max(allocations$multiplicity_of_infection),
    sum(allocations$multiplicity_of_infection == 1),
    sum(allocations$multiplicity_of_infection > 1)
  )
)

# Display summary table
knitr::kable(summary_stats, 
             col.names = c("Metric", "Value"),
             caption = "Summary Statistics from SNP-Slice Analysis")
```

| Metric              |  Value |
|:--------------------|-------:|
| Total Hosts         | 200.00 |
| Total SNPs          |  96.00 |
| Total Strains       |  52.00 |
| Mean MOI            |   2.61 |
| Max MOI             |  11.00 |
| Single Infections   | 105.00 |
| Multiple Infections |  95.00 |

Summary Statistics from SNP-Slice Analysis {.table}

## Next steps

This vignette covered the basics of using SNP-Slice with the negative
binomial model. For more:

- **Other models**: Try Poisson, binomial, or categorical models via the
  `model` argument in
  [`snp_slice()`](https://plasmogenepi.github.io/snp.slicer/reference/snp_slice.md).
- **Diagnostics**: Use
  [`?convergence_diagnostics`](https://plasmogenepi.github.io/snp.slicer/reference/convergence_diagnostics.md)
  (R-hat and ESS) and
  [`?plot_convergence`](https://plasmogenepi.github.io/snp.slicer/reference/plot_convergence.md)
  when MCMC samples are stored.
- **Your data**: Apply SNP-Slice to your own read count or categorical
  data; see
  [`?snp_slice`](https://plasmogenepi.github.io/snp.slicer/reference/snp_slice.md)
  and
  [`?load_example_results`](https://plasmogenepi.github.io/snp.slicer/reference/load_example_results.md)
  for input format and examples.
- **Parameters**: Fine-tune `alpha`, `rho`, `threshold`, `n_burnin`, and
  `gap` for your application.

## References

- SNP-Slice Resolves Mixed Infections: Simultaneously Unveiling Strain
  Haplotypes and Linking Them to Hosts
- [BioRxiv
  preprint](https://www.biorxiv.org/content/10.1101/2023.07.29.551098v2)
