# Data preparation script for SNP-Slice example data
# This script prepares the example data for inclusion in the package

library(usethis)
library(readr)

# Read the example data files
read0_file <- "data/example_read0_no_host.txt"
read1_file <- "data/example_read1_no_host.txt"

# Read the data
example_read0 <- as.matrix(read.table(read0_file, header = TRUE, sep = "\t"))
example_read1 <- as.matrix(read.table(read1_file, header = TRUE, sep = "\t"))

# Create a combined data object for easier access
example_snp_data <- list(
  read0 = example_read0,
  read1 = example_read1
)

# Use usethis to save the data to the package
usethis::use_data(example_read0, overwrite = TRUE)
usethis::use_data(example_read1, overwrite = TRUE)
usethis::use_data(example_snp_data, overwrite = TRUE)

# Print information about the data
cat("Example data prepared successfully!\n")
cat("example_read0: ", nrow(example_read0), "hosts ×", ncol(example_read0), "SNPs\n")
cat("example_read1: ", nrow(example_read1), "hosts ×", ncol(example_read1), "SNPs\n")
cat("example_snp_data: Combined data object with both read0 and read1\n")
