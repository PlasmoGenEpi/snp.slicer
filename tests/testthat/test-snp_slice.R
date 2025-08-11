test_that("snp_slice basic functionality works", {
  # Create test data
  test_data <- list(
    read1 = matrix(c(10, 5, 15, 8), nrow = 2),
    read0 = matrix(c(90, 95, 85, 92), nrow = 2)
  )
  
  # Test basic run
  result <- snp_slice(test_data, n_mcmc = 25, verbose = FALSE)
  
  # Check result structure
  expect_s3_class(result, "snp_slice_results")
  expect_true(is.matrix(result$allocation_matrix))
  expect_true(is.matrix(result$dictionary_matrix))
  expect_true(is.list(result$model_info))
  
  # Check dimensions
  expect_equal(nrow(result$allocation_matrix), 2)  # 2 hosts
  expect_equal(ncol(result$dictionary_matrix), 2)  # 4 SNPs
})

test_that("snp_slice validates input data", {
  # Test NULL data
  expect_error(snp_slice(NULL), "Data cannot be NULL or empty")
  
  # Test empty list
  expect_error(snp_slice(list()), "Data cannot be NULL or empty")
  
  # Test invalid list structure
  expect_error(snp_slice(list(read1 = matrix(1:4, 2))), 
               "List data must contain 'read1' and 'read0' elements")
  
  # Test mismatched dimensions
  expect_error(snp_slice(list(
    read1 = matrix(1:4, 2),
    read0 = matrix(1:6, 3)
  )), "read1 and read0 must have identical dimensions")
  
  # Test negative values
  expect_error(snp_slice(list(
    read1 = matrix(c(-1, 2, 3, 4), 2),
    read0 = matrix(c(5, 6, 7, 8), 2)
  )), "Read counts cannot be negative")
})

test_that("snp_slice validates parameters", {
  test_data <- list(
    read1 = matrix(c(10, 5, 15, 8), nrow = 2),
    read0 = matrix(c(90, 95, 85, 92), nrow = 2)
  )
  
  # Test invalid alpha
  expect_error(snp_slice(test_data, alpha = -1), "alpha must be a positive number")
  expect_error(snp_slice(test_data, alpha = 0), "alpha must be a positive number")
  
  # Test invalid rho
  expect_error(snp_slice(test_data, rho = -0.1), "rho must be a number between 0 and 1")
  expect_error(snp_slice(test_data, rho = 1.1), "rho must be a number between 0 and 1")
  
  # Test invalid threshold
  expect_error(snp_slice(test_data, threshold = -0.1), "threshold must be a number between 0 and 1")
  expect_error(snp_slice(test_data, threshold = 1.1), "threshold must be a number between 0 and 1")
  
  # Test invalid n_mcmc
  expect_error(snp_slice(test_data, n_mcmc = 0), "n_mcmc must be a positive integer")
  expect_error(snp_slice(test_data, n_mcmc = 1.5), "n_mcmc must be a positive integer")
})

test_that("snp_slice works with different models", {
  test_data <- list(
    read1 = matrix(c(10, 5, 15, 8), nrow = 2),
    read0 = matrix(c(90, 95, 85, 92), nrow = 2)
  )
  
  # Test negative binomial model (default)
  result_nb <- snp_slice(test_data, model = "negative_binomial", n_mcmc = 25, verbose = FALSE)
  expect_equal(result_nb$model_info$model, "negative_binomial")
  
  # Test Poisson model
  result_pois <- snp_slice(test_data, model = "poisson", n_mcmc = 25, verbose = FALSE)
  expect_equal(result_pois$model_info$model, "poisson")
  
  # Test binomial model
  result_bin <- snp_slice(test_data, model = "binomial", n_mcmc = 25, verbose = FALSE)
  expect_equal(result_bin$model_info$model, "binomial")
})

test_that("snp_slice works with categorical data", {
  # Create categorical test data
  cat_data <- matrix(c(0, 1, 0.5, 1, 0, 0.5), nrow = 2, ncol = 3)
  
  # Test categorical model
  result <- snp_slice(cat_data, model = "categorical", n_mcmc = 25, verbose = FALSE)
  expect_equal(result$model_info$model, "categorical")
  expect_equal(result$model_info$data_type, "categorical")
})

test_that("snp_slice model-specific constructors work", {
  test_data <- list(
    read1 = matrix(c(10, 5, 15, 8), nrow = 2),
    read0 = matrix(c(90, 95, 85, 92), nrow = 2)
  )
  
  # Create categorical test data
  cat_data <- matrix(c(0, 1, 0.5, 1, 0, 0.5), nrow = 2, ncol = 3)
  
  # Test categorical constructor
  result_cat <- snp_slice_categorical(cat_data, n_mcmc = 25, verbose = FALSE)
  expect_equal(result_cat$model_info$model, "categorical")
  
  # Test Poisson constructor
  result_pois <- snp_slice_poisson(test_data, n_mcmc = 25, verbose = FALSE)
  expect_equal(result_pois$model_info$model, "poisson")
  
  # Test binomial constructor
  result_bin <- snp_slice_binomial(test_data, n_mcmc = 25, verbose = FALSE)
  expect_equal(result_bin$model_info$model, "binomial")
  
  # Test negative binomial constructor
  result_nb <- snp_slice_negative_binomial(test_data, n_mcmc = 25, verbose = FALSE)
  expect_equal(result_nb$model_info$model, "negative_binomial")
})

test_that("snp_slice handles MCMC parameters correctly", {
  test_data <- list(
    read1 = matrix(c(10, 5, 15, 8), nrow = 2),
    read0 = matrix(c(90, 95, 85, 92), nrow = 2)
  )
  
  # Test with custom burnin
  result <- snp_slice(test_data, n_mcmc = 25, burnin = 5, verbose = FALSE)
  expect_equal(result$parameters$burnin, 5)
  
  # Test with gap parameter
  result <- snp_slice(test_data, n_mcmc = 25, gap = 5, verbose = FALSE)
  expect_equal(result$parameters$gap, 5)
  
  # Test with seed
  result1 <- snp_slice(test_data, n_mcmc = 25, seed = 123, verbose = FALSE)
  result2 <- snp_slice(test_data, n_mcmc = 25, seed = 123, verbose = FALSE)
  expect_equal(result1$allocation_matrix, result2$allocation_matrix)
})

test_that("snp_slice stores MCMC samples when requested", {
  test_data <- list(
    read1 = matrix(c(10, 5, 15, 8), nrow = 2),
    read0 = matrix(c(90, 95, 85, 92), nrow = 2)
  )
  
  # Test without storing samples
  result_no_samples <- snp_slice(test_data, n_mcmc = 25, store_mcmc = FALSE, verbose = FALSE)
  expect_null(result_no_samples$mcmc_samples)
  
  # Test with storing samples
  result_with_samples <- snp_slice(test_data, n_mcmc = 25, store_mcmc = TRUE, verbose = FALSE)
  expect_true(is.list(result_with_samples$mcmc_samples))
  expect_equal(length(result_with_samples$mcmc_samples), 25)
})
