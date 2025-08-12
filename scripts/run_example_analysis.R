#!/usr/bin/env Rscript

# Script to run SNP-Slice analysis on example data and save results
# This script runs once and saves the results for use in the vignette

library(snp.slicer)

# Load the example data directly from files
read0_file <- "data/example_read0_no_host.txt"
read1_file <- "data/example_read1_no_host.txt"

# Read the data
read0 <- as.matrix(read.table(read0_file, header = TRUE, sep = "\t"))
read1 <- as.matrix(read.table(read1_file, header = TRUE, sep = "\t"))

# Create data list
data <- list(read0 = read0, read1 = read1)

# Display data information
cat("Data loaded successfully!\n")
cat("Number of hosts:", nrow(read0), "\n")
cat("Number of SNPs:", ncol(read0), "\n")
cat("Read0 range:", range(read0), "\n")
cat("Read1 range:", range(read1), "\n")

# Run SNP-Slice with negative binomial model
cat("\nRunning SNP-Slice with negative binomial model...\n")
cat("This may take a few minutes...\n")

set.seed(123)  # For reproducibility
result <- snp_slice(data, 
                   model = "negative_binomial",
                   n_mcmc = 2000,
                   store_mcmc = TRUE,
                   verbose = TRUE)

# Save the results to inst/extdata/
# This is the proper location for generated data files in R packages
output_file <- "inst/extdata/example_analysis_results.rds"

# Create the directory if it doesn't exist
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

saveRDS(result, file = output_file)

cat("\nAnalysis completed and saved to:", output_file, "\n")

# Print summary information
cat("\nSummary of results:\n")
print(result)

# Extract and display key information
strains <- extract_strains(result)
allocations <- extract_allocations(result)

cat("\nKey results:\n")
cat("Number of strains identified:", strains$n_strains, "\n")
cat("Number of SNPs:", strains$n_snps, "\n")
cat("Number of hosts:", allocations$n_hosts, "\n")
cat("Multiplicity of infection (MOI) summary:\n")
print(summary(allocations$multiplicity_of_infection))

cat("\nScript completed successfully!\n")
