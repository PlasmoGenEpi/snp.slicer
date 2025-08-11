test_that("read count data files can be loaded and processed", {
  # Skip if the example files don't exist
  read1_file <- system.file("data", "example_read1_no_host.txt", package = "snp.slice")
  read0_file <- system.file("data", "example_read0_no_host.txt", package = "snp.slice")
  skip_if_not(file.exists(read1_file), "Read1 file not found")
  skip_if_not(file.exists(read0_file), "Read0 file not found")
  
  # Test that we can read the read count data files directly
  expect_no_error({
    read1_data <- utils::read.delim(read1_file, stringsAsFactors = FALSE)
    read0_data <- utils::read.delim(read0_file, stringsAsFactors = FALSE)
  })
  
  # Check dimensions
  expect_equal(nrow(read1_data), 200)  # 200 hosts
  expect_equal(ncol(read1_data), 96)   # 96 SNP sites
  expect_equal(nrow(read0_data), 200)  # 200 hosts
  expect_equal(ncol(read0_data), 96)   # 96 SNP sites
  
  # Check that data contains numeric values
  expect_true(is.numeric(as.matrix(read1_data)))
  expect_true(is.numeric(as.matrix(read0_data)))
  
  # Check that read1 values are non-negative
  expect_true(all(as.matrix(read1_data) >= 0))
  expect_true(all(as.matrix(read0_data) >= 0))
})

test_that("Poisson model works with read count data", {
  # Skip if the example files don't exist
  read1_file <- system.file("data", "example_read1_no_host.txt", package = "snp.slice")
  read0_file <- system.file("data", "example_read0_no_host.txt", package = "snp.slice")
  skip_if_not(file.exists(read1_file), "Read1 file not found")
  skip_if_not(file.exists(read0_file), "Read0 file not found")
  
  # Load read count data directly
  read1_data <- utils::read.delim(read1_file, stringsAsFactors = FALSE)
  read0_data <- utils::read.delim(read0_file, stringsAsFactors = FALSE)
  
  # Use only first 25 hosts and 25 SNPs for fast testing
  small_read1 <- as.matrix(read1_data[1:25, 1:25])
  small_read0 <- as.matrix(read0_data[1:25, 1:25])
  
  # Create data list for read count models
  data <- list(
    read1 = small_read1,
    read0 = small_read0
  )
  
  # Test running SNP-Slice with the Poisson model
  expect_no_error({
    result <- snp_slice_poisson(
        data,
        n_mcmc = 25,  # Short run for testing
        verbose = FALSE,
        store_mcmc = TRUE
    )
  })
  
  # Check that the result has the expected structure
  expect_s3_class(result, "snp_slice_results")
  expect_equal(result$model_info$model, "poisson")
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
  expect_equal(allocations$n_hosts, nrow(small_read1))  # Same number of hosts
})

test_that("Negative Binomial model works with read count data", {
  # Skip if the example files don't exist
  read1_file <- system.file("data", "example_read1_no_host.txt", package = "snp.slice")
  read0_file <- system.file("data", "example_read0_no_host.txt", package = "snp.slice")
  skip_if_not(file.exists(read1_file), "Read1 file not found")
  skip_if_not(file.exists(read0_file), "Read0 file not found")
  
  # Load read count data directly
  read1_data <- utils::read.delim(read1_file, stringsAsFactors = FALSE)
  read0_data <- utils::read.delim(read0_file, stringsAsFactors = FALSE)
  
  # Use only first 25 hosts and 25 SNPs for fast testing
  small_read1 <- as.matrix(read1_data[1:25, 1:25])
  small_read0 <- as.matrix(read0_data[1:25, 1:25])
  
  # Create data list for read count models
  data <- list(
    read1 = small_read1,
    read0 = small_read0
  )
  
  # Test running SNP-Slice with the Negative Binomial model
  expect_no_error({
    result <- snp_slice_negative_binomial(list(
      read1 = small_read1,
      read0 = small_read0
    ), 
                                         n_mcmc = 25,  # Short run for testing
                                         verbose = FALSE,
                                         store_mcmc = TRUE)
  })
  
  # Check that the result has the expected structure
  expect_s3_class(result, "snp_slice_results")
  expect_equal(result$model_info$model, "negative_binomial")
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
  expect_equal(allocations$n_hosts, nrow(small_read1))  # Same number of hosts
})

test_that("Binomial model works with read count models", {
  # Skip if the example files don't exist
  read1_file <- system.file("data", "example_read1_no_host.txt", package = "snp.slice")
  read0_file <- system.file("data", "example_read0_no_host.txt", package = "snp.slice")
  skip_if_not(file.exists(read1_file), "Read1 file not found")
  skip_if_not(file.exists(read0_file), "Read0 file not found")
  
  # Load read count data directly
  read1_data <- utils::read.delim(read1_file, stringsAsFactors = FALSE)
  read0_data <- utils::read.delim(read0_file, stringsAsFactors = FALSE)
  
  # Use only first 25 hosts and 25 SNPs for fast testing
  small_read1 <- as.matrix(read1_data[1:25, 1:25])
  small_read0 <- as.matrix(read0_data[1:25, 1:25])
  
  # Create data list for read count models
  data <- list(
    read1 = small_read1,
    read0 = small_read0
  )
  
  # Test running SNP-Slice with the Binomial model
  expect_no_error({
    result <- snp_slice_binomial(list(
      read1 = small_read1,
      read0 = small_read0
    ), 
                                n_mcmc = 25,  # Short run for testing
                                verbose = FALSE,
                                store_mcmc = TRUE)
  })
  
  # Check that the result has the expected structure
  expect_s3_class(result, "snp_slice_results")
  expect_equal(result$model_info$model, "binomial")
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
  expect_equal(allocations$n_hosts, nrow(small_read1))  # Same number of hosts
})

test_that("read count data preprocessing works correctly", {
  # Skip if the example files don't exist
  read1_file <- system.file("data", "example_read1_no_host.txt", package = "snp.slice")
  read0_file <- system.file("data", "example_read0_no_host.txt", package = "snp.slice")
  skip_if_not(file.exists(read1_file), "Read1 file not found")
  skip_if_not(file.exists(read0_file), "Read0 file not found")
  
  # Load read count data directly
  read1_data <- utils::read.delim(read1_file, stringsAsFactors = FALSE)
  read0_data <- utils::read.delim(read0_file, stringsAsFactors = FALSE)
  
  # Use only first 25 hosts and 25 SNPs for fast testing
  small_read1 <- as.matrix(read1_data[1:25, 1:25])
  small_read0 <- as.matrix(read0_data[1:25, 1:25])
  
  # Create data list for read count models
  data <- list(
    read1 = small_read1,
    read0 = small_read0
  )
  
  # Test preprocessing for each model
  models <- c("poisson", "negative_binomial", "binomial")
  
  for (model in models) {
    expect_no_error({
      processed_data <- preprocess_data(list(
        read1 = small_read1,
        read0 = small_read0
      ), model)
    })
    
    # Check processed data structure
    expect_equal(processed_data$data_type, "read_counts")
    expect_equal(processed_data$model, model)
    expect_equal(processed_data$N, 25)  # 25 hosts
    expect_equal(processed_data$P, 25)  # 25 SNP sites
    
    # Check that y and r matrices have correct dimensions
    expect_equal(nrow(processed_data$y), 25)
    expect_equal(ncol(processed_data$y), 25)
    expect_equal(nrow(processed_data$r), 25)
    expect_equal(ncol(processed_data$r), 25)
    
    # Check that r = read1 + read0
    expected_r <- small_read1 + small_read0
    expect_equal(processed_data$r, expected_r)
    
    # Check that y = read1
    expect_equal(processed_data$y, small_read1)
  }
})

test_that("read count models converge and produce reasonable results", {
  # Skip if the example files don't exist
  read1_file <- system.file("data", "example_read1_no_host.txt", package = "snp.slice")
  read0_file <- system.file("data", "example_read0_no_host.txt", package = "snp.slice")
  skip_if_not(file.exists(read1_file), "Read1 file not found")
  skip_if_not(file.exists(read0_file), "Read0 file not found")
  
  # Load read count data directly
  read1_data <- utils::read.delim(read1_file, stringsAsFactors = FALSE)
  read0_data <- utils::read.delim(read0_file, stringsAsFactors = FALSE)
  
  # Use only first 25 hosts and 25 SNPs for fast testing
  small_read1 <- as.matrix(read1_data[1:25, 1:25])
  small_read0 <- as.matrix(read0_data[1:25, 1:25])
  
  # Create data list for read count models
  data <- list(
    read1 = small_read1,
    read0 = small_read0
  )
  
  # Test each model with a longer run to check convergence
  models <- list(
    poisson = snp_slice_poisson,
    negative_binomial = snp_slice_negative_binomial,
    binomial = snp_slice_binomial
  )
  
  for (model_name in names(models)) {
    model_func <- models[[model_name]]
    
    expect_no_error({
      result <- model_func(list(
        read1 = small_read1,
        read0 = small_read0
      ), 
                          n_mcmc = 25,  # Short run for convergence testing
                          verbose = FALSE,
                          store_mcmc = TRUE)
    })
    
    # Check that the model converged (log-posterior is finite)
    expect_true(is.finite(result$diagnostics$final_logpost))
    expect_true(is.finite(result$diagnostics$map_logpost))
    
    # Check that we have some active strains
    expect_true(result$diagnostics$final_kstar > 0)
    
    # Check that the number of active strains is reasonable
    # (should be less than the number of hosts)
    expect_true(result$diagnostics$final_kstar <= nrow(small_read1))
    
    # Check that we have some active strains
    expect_true(result$diagnostics$final_kstar > 0)
    expect_true(result$diagnostics$final_kstar <= nrow(small_read1))
  }
})

test_that("read count models handle different parameter settings", {
  # Skip if the example files don't exist
  read1_file <- system.file("data", "example_read1_no_host.txt", package = "snp.slice")
  read0_file <- system.file("data", "example_read0_no_host.txt", package = "snp.slice")
  skip_if_not(file.exists(read1_file), "Read1 file not found")
  skip_if_not(file.exists(read0_file), "Read0 file not found")
  
  # Load read count data directly
  read1_data <- utils::read.delim(read1_file, stringsAsFactors = FALSE)
  read0_data <- utils::read.delim(read0_file, stringsAsFactors = FALSE)
  
  # Use only first 25 hosts and 25 SNPs for fast testing
  small_read1 <- as.matrix(read1_data[1:25, 1:25])
  small_read0 <- as.matrix(read0_data[1:25, 1:25])
  
  # Create data list for read count models
  data <- list(
    read1 = small_read1,
    read0 = small_read0
  )
  
  # Test with different alpha and rho values
  test_params <- list(
    list(alpha = 1.0, rho = 0.3),
    list(alpha = 5.0, rho = 0.7),
    list(alpha = 2.6, rho = 0.5)  # Default values
  )
  
  for (params in test_params) {
    expect_no_error({
      result <- snp_slice_poisson(list(
        read1 = small_read1,
        read0 = small_read0
      ), 
                                 n_mcmc = 25,  # Short run for testing
                                 alpha = params$alpha,
                                 rho = params$rho,
                                 verbose = FALSE)
    })
    
    # Check that the model ran successfully
    expect_true(is.finite(result$diagnostics$final_logpost))
    expect_true(result$diagnostics$final_kstar > 0)
  }
})
