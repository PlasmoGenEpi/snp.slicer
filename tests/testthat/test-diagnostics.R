test_that("extract_strains works correctly", {
  # Create test result
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3)
  )
  class(test_result) <- "snp_slice_results"
  
  # Test extraction
  strains <- extract_strains(test_result)
  
  # Check structure
  expect_true(is.list(strains))
  expect_equal(strains$n_strains, 2)
  expect_equal(strains$n_snps, 3)
  expect_equal(length(strains$strain_names), 2)
  expect_true(is.matrix(strains$dictionary))
  expect_equal(dim(strains$dictionary), c(2, 3))
})

test_that("extract_allocations works correctly", {
  # Create test result
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3)
  )
  class(test_result) <- "snp_slice_results"
  
  # Test extraction
  allocations <- extract_allocations(test_result)
  
  # Check structure
  expect_true(is.list(allocations))
  expect_equal(allocations$n_hosts, 2)
  expect_equal(allocations$n_strains, 2)
  expect_equal(length(allocations$host_names), 2)
  expect_equal(length(allocations$strain_names), 2)
  expect_true(is.matrix(allocations$allocation_matrix))
  expect_equal(dim(allocations$allocation_matrix), c(2, 2))
  
  # Check MOI calculation
  expect_equal(allocations$multiplicity_of_infection, c(2, 1))
})

test_that("extract functions validate input", {
  # Test with non-snp_slice_results object
  expect_error(extract_strains(list()), "Input must be an snp_slice_results object")
  expect_error(extract_allocations(list()), "Input must be an snp_slice_results object")
})

test_that("summary.snp_slice_results works", {
  # Create test result
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3, data_type = "read_counts"),
    convergence = list(gap_converged = TRUE, iterations_run = 100),
    diagnostics = list(
      final_logpost = -10.5,
      map_logpost = -9.8,
      final_kstar = 2,
      map_kstar = 2
    )
  )
  class(test_result) <- "snp_slice_results"
  
  # Test summary (capture output)
  output <- capture.output(summary(test_result))
  
  # Check that summary produces output
  expect_true(length(output) > 0)
  expect_true(any(grepl("SNP-Slice Results Summary", output)))
  expect_true(any(grepl("negative_binomial", output)))
  expect_true(any(grepl("2 hosts", output)))
  expect_true(any(grepl("3 SNPs", output)))
})

test_that("print.snp_slice_results works", {
  # Create test result
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3),
    convergence = list(gap_converged = TRUE, iterations_run = 100)
  )
  class(test_result) <- "snp_slice_results"
  
  # Test print (capture output)
  output <- capture.output(print(test_result))
  
  # Check that print produces output
  expect_true(length(output) > 0)
  expect_true(any(grepl("SNP-Slice Results", output)))
  expect_true(any(grepl("negative_binomial", output)))
  expect_true(any(grepl("2 hosts", output)))
#   expect_true(any(grepl("2 strains", output)))
})

test_that("plot_convergence validates input", {
  # Test with non-snp_slice_results object
  expect_error(plot_convergence(list()), "Input must be an snp_slice_results object")
  
  # Test with result without MCMC samples
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3)
  )
  class(test_result) <- "snp_slice_results"
  
  expect_error(plot_convergence(test_result), "MCMC samples not stored")
  
  # Test with invalid plot type
  test_result$mcmc_samples <- list(
    list(logpost = -10, kstar = 2, A = matrix(1, 2, 2)),
    list(logpost = -9, kstar = 2, A = matrix(1, 2, 2))
  )
  
  expect_error(plot_convergence(test_result, type = "invalid"), 
               "Invalid plot type")
})

test_that("effective_sample_size validates input", {
  # Test with non-snp_slice_results object
  expect_error(effective_sample_size(list()), "Input must be an snp_slice_results object")
  
  # Test with result without MCMC samples
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3)
  )
  class(test_result) <- "snp_slice_results"
  
  expect_error(effective_sample_size(test_result), "MCMC samples not stored")
  
  # Test with invalid parameter
  test_result$mcmc_samples <- list(
    list(logpost = -10, kstar = 2, A = matrix(1, 2, 2)),
    list(logpost = -9, kstar = 2, A = matrix(1, 2, 2))
  )
  
  expect_error(effective_sample_size(test_result, parameter = "invalid"), 
               "Invalid parameter")
})

test_that("effective_sample_size works with scalar parameters", {
  # Create test result with MCMC samples
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3),
    mcmc_samples = list(
      list(logpost = -10, kstar = 2, ktrunc = 3, A = matrix(1, 2, 2)),
      list(logpost = -9, kstar = 2, ktrunc = 3, A = matrix(1, 2, 2)),
      list(logpost = -8, kstar = 3, ktrunc = 4, A = matrix(1, 2, 2))
    )
  )
  class(test_result) <- "snp_slice_results"
  
  # Test logpost ESS
  ess_logpost <- effective_sample_size(test_result, parameter = "logpost")
  expect_true(is.list(ess_logpost))
  expect_equal(ess_logpost$parameter, "logpost")
  expect_equal(ess_logpost$method, "autocorrelation")
  expect_equal(ess_logpost$n_samples, 3)
  expect_true(ess_logpost$ess > 0)
  expect_true(ess_logpost$ess <= 3)
  expect_true("efficiency" %in% names(ess_logpost))
  
  # Test kstar ESS
  ess_kstar <- effective_sample_size(test_result, parameter = "kstar")
  expect_equal(ess_kstar$parameter, "kstar")
  expect_true(ess_kstar$ess > 0)
  
  # Test ktrunc ESS
  ess_ktrunc <- effective_sample_size(test_result, parameter = "ktrunc")
  expect_equal(ess_ktrunc$parameter, "ktrunc")
  expect_true(ess_ktrunc$ess > 0)
  
  # Test n_strains ESS
  ess_n_strains <- effective_sample_size(test_result, parameter = "n_strains")
  expect_equal(ess_n_strains$parameter, "n_strains")
  expect_true(ess_n_strains$ess > 0)
})

test_that("effective_sample_size works with different methods", {
  # Create test result with MCMC samples
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3),
    mcmc_samples = list(
      list(logpost = -10, kstar = 2, A = matrix(1, 2, 2)),
      list(logpost = -9, kstar = 2, A = matrix(1, 2, 2)),
      list(logpost = -8, kstar = 3, A = matrix(1, 2, 2))
    )
  )
  class(test_result) <- "snp_slice_results"
  
  # Test autocorrelation method (default)
  ess_auto <- effective_sample_size(test_result, parameter = "logpost", method = "autocorrelation")
  expect_equal(ess_auto$method, "autocorrelation")
  
  # Test batch means method
  ess_batch <- effective_sample_size(test_result, parameter = "logpost", method = "batch_means")
  expect_equal(ess_batch$method, "batch_means")
  
  # Test spectral method
  ess_spec <- effective_sample_size(test_result, parameter = "logpost", method = "spectral")
  expect_equal(ess_spec$method, "spectral")
})

test_that("effective_sample_size works with matrix parameters", {
  # Create test result with MCMC samples including A and D matrices
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3),
    mcmc_samples = list(
      list(logpost = -10, kstar = 2, A = matrix(c(1, 0, 1, 1), nrow = 2), 
           D = matrix(c(1, 0, 0, 1), nrow = 2)),
      list(logpost = -9, kstar = 2, A = matrix(c(1, 1, 0, 1), nrow = 2), 
           D = matrix(c(1, 1, 0, 0), nrow = 2)),
      list(logpost = -8, kstar = 3, A = matrix(c(0, 1, 1, 0), nrow = 2), 
           D = matrix(c(0, 1, 1, 1), nrow = 2))
    )
  )
  class(test_result) <- "snp_slice_results"
  
  # Test A matrix ESS
  ess_A <- effective_sample_size(test_result, parameter = "A")
  expect_equal(ess_A$parameter, "A")
  expect_true("overall_ess" %in% names(ess_A))
  expect_true("matrix_dimensions" %in% names(ess_A))
  expect_true("components" %in% names(ess_A))
  expect_equal(ess_A$matrix_dimensions, c(2, 2))
  expect_true(ess_A$overall_ess > 0)
  
  # Test D matrix ESS
  ess_D <- effective_sample_size(test_result, parameter = "D")
  expect_equal(ess_D$parameter, "D")
  expect_true("overall_ess" %in% names(ess_D))
  expect_true("matrix_dimensions" %in% names(ess_D))
  expect_equal(ess_D$matrix_dimensions, c(2, 2))
  expect_true(ess_D$overall_ess > 0)
})

test_that("effective_sample_size works with mu parameter", {
  # Create test result with MCMC samples including mu
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3),
    mcmc_samples = list(
      list(logpost = -10, kstar = 2, mu = c(0.8, 0.6, 0.4)),
      list(logpost = -9, kstar = 2, mu = c(0.7, 0.5, 0.3)),
      list(logpost = -8, kstar = 3, mu = c(0.9, 0.7, 0.5, 0.2))
    )
  )
  class(test_result) <- "snp_slice_results"
  
  # Test mu ESS
  ess_mu <- effective_sample_size(test_result, parameter = "mu")
  expect_equal(ess_mu$parameter, "mu")
  expect_true("components" %in% names(ess_mu))
  expect_true(length(ess_mu$components) > 0)
  
  # Check that each component has ESS
  for (comp_name in names(ess_mu$components)) {
    comp <- ess_mu$components[[comp_name]]
    expect_true(comp$ess > 0)
    expect_true("n_samples" %in% names(comp))
  }
})

test_that("effective_sample_size works with 'all' parameter", {
  # Create test result with comprehensive MCMC samples
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3),
    mcmc_samples = list(
      list(logpost = -10, kstar = 2, ktrunc = 3, 
           A = matrix(c(1, 0, 1, 1), nrow = 2), 
           D = matrix(c(1, 0, 0, 1), nrow = 2),
           mu = c(0.8, 0.6, 0.4)),
      list(logpost = -9, kstar = 2, ktrunc = 3, 
           A = matrix(c(1, 1, 0, 1), nrow = 2), 
           D = matrix(c(1, 1, 0, 0), nrow = 2),
           mu = c(0.7, 0.5, 0.3)),
      list(logpost = -8, kstar = 3, ktrunc = 4, 
           A = matrix(c(0, 1, 1, 0), nrow = 2), 
           D = matrix(c(0, 1, 1, 1), nrow = 2),
           mu = c(0.9, 0.7, 0.5, 0.2))
    )
  )
  class(test_result) <- "snp_slice_results"
  
  # Test 'all' parameter
  ess_all <- effective_sample_size(test_result, parameter = "all")
  expect_true(is.list(ess_all))
  expect_true(length(ess_all) > 0)
  
  # Check that expected parameters are present
  expected_params <- c("logpost", "kstar", "ktrunc", "mu", "A", "D")
  for (param in expected_params) {
    expect_true(param %in% names(ess_all))
  }
  
  # Check that each parameter has valid ESS results
  for (param in names(ess_all)) {
    param_result <- ess_all[[param]]
    expect_true(is.list(param_result))
    expect_equal(param_result$parameter, param)
  }
})

test_that("effective_sample_size handles edge cases", {
  # Test with single sample
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(c(1, 0, 0, 1, 1, 0), nrow = 2, ncol = 3),
    model_info = list(model = "negative_binomial", N = 2, P = 3),
    mcmc_samples = list(
      list(logpost = -10, kstar = 2, A = matrix(1, 2, 2))
    )
  )
  class(test_result) <- "snp_slice_results"
  
  # Should handle single sample gracefully
  expect_warning(ess_single <- effective_sample_size(test_result, parameter = "logpost"), 
                "Insufficient samples for ESS calculation")
  expect_equal(ess_single$ess, 1)
  expect_equal(ess_single$n_samples, 1)
  
  # Test with constant values (should still work)
  test_result$mcmc_samples <- list(
    list(logpost = -10, kstar = 2, A = matrix(1, 2, 2)),
    list(logpost = -10, kstar = 2, A = matrix(1, 2, 2)),
    list(logpost = -10, kstar = 2, A = matrix(1, 2, 2))
  )
  
  ess_constant <- effective_sample_size(test_result, parameter = "logpost")
  expect_true(ess_constant$ess > 0)
  expect_true(ess_constant$ess <= 3)
})
