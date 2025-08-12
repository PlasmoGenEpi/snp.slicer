test_that("example categorical file can be processed and run", {
  # Skip if the example file doesn't exist
  example_file <- system.file("data", "example_cat.txt", package = "snp.slicer")
  skip_if_not(file.exists(example_file), "Example file not found")
  
  # Test that we can read the file
  expect_no_error({
    data <- read_snp_data(example_file)
  })
  
  # Check the structure of the loaded data
  expect_true(is.matrix(data))
  expect_equal(ncol(data), 96)  # Should have 96 SNP sites (removing host_id and strain_id)
  expect_true(all(data %in% c(0, 1)))  # Should be binary categorical data
  
  # Use only first 25 hosts and 25 SNPs for fast testing
  small_data <- data[1:25, 1:25]
  
  # Test running SNP-Slice with the categorical model
  expect_no_error({
    result <- snp_slice_categorical(small_data, 
                                   n_mcmc = 25,  # Short run for testing
                                   verbose = FALSE,
                                   store_mcmc = TRUE)
  })
  
  # Check that the result has the expected structure
  expect_s3_class(result, "snp_slice_results")
  expect_equal(result$model_info$model, "categorical")
  expect_true(length(result$mcmc_samples) > 0)
  expect_true("A" %in% names(result$mcmc_samples[[1]]))
  expect_true("D" %in% names(result$mcmc_samples[[1]]))
  expect_true("mu" %in% names(result$mcmc_samples[[1]]))
  
  # Test that we can extract strains and allocations
  expect_no_error({
    strains <- extract_strains(result)
    allocations <- extract_allocations(result)
  })
  
  # Check that extracted data has reasonable dimensions
  expect_true(is.list(strains))
  expect_true(is.list(allocations))
  expect_equal(allocations$n_hosts, nrow(small_data))  # Same number of hosts
})

test_that("example file data preprocessing works correctly", {
  # Skip if the example file doesn't exist
  example_file <- system.file("data", "example_cat.txt", package = "snp.slicer")
  skip_if_not(file.exists(example_file), "Example file not found")
  
  # Test the preprocessing function directly
  expect_no_error({
    processed_data <- preprocess_data(example_file, model = "categorical")
  })
  
  # Check that preprocessing removed the ID columns
  expect_equal(ncol(processed_data$y), 96)  # Should have 96 SNP sites
  expect_equal(nrow(processed_data$y), nrow(read.delim(example_file)))  # Same number of hosts
  
  # Load data for small subset testing
  data <- read_snp_data(example_file)
  small_data <- data[1:25, 1:25]
  
  # Test with small subset for faster testing
  small_processed_data <- preprocess_data(small_data, model = "categorical")
  expect_equal(ncol(small_processed_data$y), 25)  # Should have 25 SNP sites
  expect_equal(nrow(small_processed_data$y), 25)  # Should have 25 hosts
  
  # Check that the data is properly formatted for categorical model
  expect_true(all(processed_data$y %in% c(0, 1)))
  expect_equal(processed_data$model, "categorical")
})
