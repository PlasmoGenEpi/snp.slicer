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

test_that("calculate_individual_coi with use_map = TRUE (MAP only)", {
  A <- matrix(c(1, 1, 1, 0, 0, 1), nrow = 3, ncol = 2, byrow = TRUE)
  test_result <- list(
    allocation_matrix = A,
    dictionary_matrix = matrix(0, nrow = 2, ncol = 2),
    model_info = list(processed_data = list(specimen_ids = NULL))
  )
  class(test_result) <- "snp_slice_results"

  coi <- calculate_individual_coi(test_result, use_map = TRUE)

  expect_equal(nrow(coi), 3)
  expect_equal(coi$host_index, 1:3)
  expect_equal(coi$coi_estimate, c(2, 1, 1))
  expect_true(all(is.na(coi$coi_sd)))
  expect_true(all(is.na(coi$coi_lower)))
  expect_true(all(is.na(coi$coi_upper)))
})

test_that("calculate_individual_coi with MCMC samples (use_map = FALSE)", {
  A <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2)
  test_result <- list(
    allocation_matrix = A,
    dictionary_matrix = matrix(0, nrow = 2, ncol = 2),
    model_info = list(processed_data = list(specimen_ids = NULL)),
    mcmc_samples = list(
      list(A = A), list(A = A), list(A = A), list(A = A), list(A = A)
    )
  )
  class(test_result) <- "snp_slice_results"

  coi <- calculate_individual_coi(test_result, use_map = FALSE, n_samples = 4, interval = 0.95)

  expect_equal(nrow(coi), 2)
  expect_equal(coi$coi_estimate, c(2, 1))
  expect_equal(coi$coi_sd, c(0, 0))
  expect_equal(coi$coi_lower, c(2, 1))
  expect_equal(coi$coi_upper, c(2, 1))
})

test_that("calculate_individual_coi errors when use_map = FALSE and no MCMC samples", {
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(0, nrow = 2, ncol = 2),
    model_info = list(processed_data = list(specimen_ids = NULL))
  )
  class(test_result) <- "snp_slice_results"

  expect_error(
    calculate_individual_coi(test_result, use_map = FALSE),
    "MCMC samples not available"
  )
})

test_that("calculate_individual_coi validates input", {
  test_result <- list(
    allocation_matrix = matrix(c(1, 0, 1, 1), nrow = 2),
    dictionary_matrix = matrix(0, nrow = 2, ncol = 2),
    model_info = list(processed_data = list(specimen_ids = NULL))
  )
  class(test_result) <- "snp_slice_results"

  expect_error(calculate_individual_coi(list()), "results must be an snp_slice_results object")
  expect_error(calculate_individual_coi(test_result, use_map = "yes"), "use_map must be a single logical")
  expect_error(calculate_individual_coi(test_result, n_samples = 0), "n_samples must be a positive integer")
  expect_error(calculate_individual_coi(test_result, interval = 1), "interval must be a number between 0 and 1")
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
    list(logpost = -10, kstar = 2, A = matrix(1, 2, 2), iteration = 1),
    list(logpost = -9, kstar = 2, A = matrix(1, 2, 2), iteration = 2)
  )
  
  expect_error(plot_convergence(test_result, type = "invalid"), 
               "Invalid plot type")
})
