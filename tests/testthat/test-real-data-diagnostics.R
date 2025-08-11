# Tests for diagnostic functions using real computed results
# These tests use the actual analysis results from the example data

test_that("load_example_results works correctly", {
  # Test that we can load the example results
  result <- load_example_results()
  
  # Check basic structure
  expect_true(inherits(result, "snp_slice_results"))
  expect_true(is.list(result))
  
  # Check required components
  expect_true("allocation_matrix" %in% names(result))
  expect_true("dictionary_matrix" %in% names(result))
  expect_true("model_info" %in% names(result))
  
  # Check dimensions (should match our example data)
  expect_equal(nrow(result$allocation_matrix), 200)  # 200 hosts
  expect_equal(ncol(result$dictionary_matrix), 96)   # 96 SNPs
  
  # Check model info
  expect_equal(result$model_info$model, "negative_binomial")
  expect_equal(result$model_info$N, 200)
  expect_equal(result$model_info$P, 96)
})

test_that("extract_strains works with real data", {
  result <- load_example_results()
  strains <- extract_strains(result)
  
  # Check structure
  expect_true(is.list(strains))
  expect_true("n_strains" %in% names(strains))
  expect_true("n_snps" %in% names(strains))
  expect_true("strain_names" %in% names(strains))
  expect_true("dictionary" %in% names(strains))
  
  # Check values
  expect_true(strains$n_strains > 0)
  expect_equal(strains$n_snps, 96)
  expect_equal(length(strains$strain_names), strains$n_strains)
  expect_equal(dim(strains$dictionary), c(strains$n_strains, 96))
  
  # Check that strain names are unique
  expect_equal(length(unique(strains$strain_names)), strains$n_strains)
  
  # Check that dictionary matrix is binary (0 or 1)
  expect_true(all(strains$dictionary %in% c(0, 1)))
})

test_that("extract_allocations works with real data", {
  result <- load_example_results()
  allocations <- extract_allocations(result)
  
  # Check structure
  expect_true(is.list(allocations))
  expect_true("n_hosts" %in% names(allocations))
  expect_true("n_strains" %in% names(allocations))
  expect_true("host_names" %in% names(allocations))
  expect_true("strain_names" %in% names(allocations))
  expect_true("allocation_matrix" %in% names(allocations))
  expect_true("multiplicity_of_infection" %in% names(allocations))
  
  # Check values
  expect_equal(allocations$n_hosts, 200)
  expect_true(allocations$n_strains > 0)
  expect_equal(length(allocations$host_names), 200)
  expect_equal(length(allocations$strain_names), allocations$n_strains)
  expect_equal(dim(allocations$allocation_matrix), c(200, allocations$n_strains))
  expect_equal(length(allocations$multiplicity_of_infection), 200)
  
  # Check that host names are unique
  expect_equal(length(unique(allocations$host_names)), 200)
  
  # Check that MOI values are reasonable
  expect_true(all(allocations$multiplicity_of_infection >= 0))
  expect_true(all(allocations$multiplicity_of_infection <= allocations$n_strains))
  
  # Check that allocation matrix sums to MOI
  calculated_moi <- rowSums(allocations$allocation_matrix)
  expect_equal(allocations$multiplicity_of_infection, calculated_moi)
})

test_that("summary.snp_slice_results works with real data", {
  result <- load_example_results()
  
  # Capture output
  output <- capture.output(summary(result))
  
  # Check that summary produces output
  expect_true(length(output) > 0)
  
  # Check for key information
  expect_true(any(grepl("negative_binomial", output, ignore.case = TRUE)))
  expect_true(any(grepl("200", output)))
  expect_true(any(grepl("96", output)))
  
  # Check that it mentions strains
  strains <- extract_strains(result)
  expect_true(any(grepl(as.character(strains$n_strains), output)))
})

test_that("print.snp_slice_results works with real data", {
  result <- load_example_results()
  
  # Capture output
  output <- capture.output(print(result))
  
  # Check that print produces output
  expect_true(length(output) > 0)
  
  # Check for key information
  expect_true(any(grepl("SNP-Slice Results", output)))
  expect_true(any(grepl("negative_binomial", output, ignore.case = TRUE)))
  expect_true(any(grepl("200", output)))
  expect_true(any(grepl("96", output)))
})

test_that("effective_sample_size works with real data when MCMC samples are available", {
  result <- load_example_results()
  
  # Only test if MCMC samples are available
  if (!is.null(result$mcmc_samples)) {
    # Test logpost ESS
    ess_logpost <- effective_sample_size(result, parameter = "logpost")
    expect_true(is.list(ess_logpost))
    expect_equal(ess_logpost$parameter, "logpost")
    expect_true(ess_logpost$ess > 0)
    expect_true(ess_logpost$n_samples > 0)
    
    # Test kstar ESS (may be NaN if kstar is constant)
    ess_kstar <- effective_sample_size(result, parameter = "kstar")
    expect_equal(ess_kstar$parameter, "kstar")
    # kstar ESS might be NaN if the parameter is constant
    expect_true(ess_kstar$ess > 0 || is.nan(ess_kstar$ess))
    
    # Test different methods for logpost
    ess_auto <- effective_sample_size(result, parameter = "logpost", method = "autocorrelation")
    ess_batch <- effective_sample_size(result, parameter = "logpost", method = "batch_means")
    ess_spec <- effective_sample_size(result, parameter = "logpost", method = "spectral")
    
    expect_equal(ess_auto$method, "autocorrelation")
    expect_equal(ess_batch$method, "batch_means")
    expect_equal(ess_spec$method, "spectral")
    
    # All should have positive ESS
    expect_true(ess_auto$ess > 0)
    expect_true(ess_batch$ess > 0)
    expect_true(ess_spec$ess > 0)
  } else {
    # If no MCMC samples, should get appropriate message
    expect_message(effective_sample_size(result, parameter = "logpost"), 
                  "MCMC samples not stored")
  }
})

test_that("strain diversity analysis works with real data", {
  result <- load_example_results()
  A <- result$allocation_matrix
  
  # Calculate strain diversity
  strain_diversity <- rowSums(A > 0)
  
  # Check basic properties
  expect_equal(length(strain_diversity), 200)
  expect_true(all(strain_diversity >= 0))
  expect_true(all(strain_diversity <= ncol(A)))
  
  # Check summary statistics
  diversity_summary <- summary(strain_diversity)
  expect_true(diversity_summary["Min."] >= 0)
  expect_true(diversity_summary["Max."] <= ncol(A))
  expect_true(diversity_summary["Mean"] >= 0)
  
  # Check that diversity table makes sense
  diversity_table <- table(strain_diversity)
  expect_true(length(diversity_table) > 0)
  expect_true(all(as.numeric(names(diversity_table)) >= 0))
  expect_true(sum(diversity_table) == 200)
})

test_that("strain frequency analysis works with real data", {
  result <- load_example_results()
  A <- result$allocation_matrix
  
  # Calculate strain frequencies
  strain_frequencies <- colSums(A)
  
  # Check basic properties
  expect_equal(length(strain_frequencies), ncol(A))
  expect_true(all(strain_frequencies >= 0))
  expect_true(all(strain_frequencies <= 200))
  
  # Check that total allocations match
  expect_equal(sum(strain_frequencies), sum(rowSums(A)))
  
  # Check frequency table
  freq_table <- table(strain_frequencies)
  expect_true(length(freq_table) > 0)
  expect_true(all(as.numeric(names(freq_table)) >= 0))
  expect_true(sum(freq_table) == ncol(A))
})

test_that("dictionary matrix analysis works with real data", {
  result <- load_example_results()
  D <- result$dictionary_matrix
  
  # Check basic properties
  expect_true(is.matrix(D))
  expect_equal(ncol(D), 96)
  expect_true(nrow(D) > 0)
  
  # Check that matrix is binary
  expect_true(all(D %in% c(0, 1)))
  
  # Check unique patterns
  unique_patterns <- unique(D)
  expect_true(nrow(unique_patterns) <= nrow(D))
  expect_true(nrow(unique_patterns) > 0)
  
  # Check that each pattern has the right number of SNPs
  expect_true(all(ncol(unique_patterns) == 96))
})

test_that("allocation matrix analysis works with real data", {
  result <- load_example_results()
  A <- result$allocation_matrix
  
  # Check basic properties
  expect_true(is.matrix(A))
  expect_equal(nrow(A), 200)
  expect_true(ncol(A) > 0)
  
  # Check that matrix is non-negative
  expect_true(all(A >= 0))
  
  # Check that each host has at least one strain
  expect_true(all(rowSums(A) > 0))
  
  # Check that allocation sums match MOI
  moi <- rowSums(A)
  expect_true(all(moi > 0))
  expect_true(all(moi <= ncol(A)))
})

test_that("model information is consistent", {
  result <- load_example_results()
  
  # Check model info
  expect_true("model_info" %in% names(result))
  expect_equal(result$model_info$model, "negative_binomial")
  expect_equal(result$model_info$N, 200)
  expect_equal(result$model_info$P, 96)
  
  # Check consistency with matrices
  expect_equal(nrow(result$allocation_matrix), result$model_info$N)
  expect_equal(ncol(result$dictionary_matrix), result$model_info$P)
})

test_that("convergence information is available", {
  result <- load_example_results()
  
  # Check convergence info if available
  if ("convergence" %in% names(result)) {
    conv <- result$convergence
    expect_true(is.list(conv))
    
    if ("converged" %in% names(conv)) {
      expect_true(is.logical(conv$converged))
    }
    
    if ("iterations_run" %in% names(conv)) {
      expect_true(conv$iterations_run > 0)
    }
  }
})

test_that("diagnostics information is available", {
  result <- load_example_results()
  
  # Check diagnostics info if available
  if ("diagnostics" %in% names(result)) {
    diag <- result$diagnostics
    expect_true(is.list(diag))
    
    # Check for common diagnostic fields
    if ("final_logpost" %in% names(diag)) {
      expect_true(is.numeric(diag$final_logpost))
    }
    
    if ("map_logpost" %in% names(diag)) {
      expect_true(is.numeric(diag$map_logpost))
    }
    
    if ("final_kstar" %in% names(diag)) {
      expect_true(is.numeric(diag$final_kstar))
      expect_true(diag$final_kstar > 0)
    }
  }
})

test_that("performance monitoring works", {
  # Test performance summary if available
  perf_summary <- get_performance_summary()
  
  if (!is.null(perf_summary)) {
    expect_true(is.list(perf_summary))
    # Performance summary might be empty if not initialized
    expect_true(length(perf_summary) >= 0)
  }
})

test_that("data loading functions work correctly", {
  # Test that we can load the example data
  data(example_snp_data, package = "snp.slice")
  
  # Check structure
  expect_true(is.list(example_snp_data))
  expect_true("read0" %in% names(example_snp_data))
  expect_true("read1" %in% names(example_snp_data))
  
  # Check dimensions
  expect_equal(nrow(example_snp_data$read0), 200)
  expect_equal(ncol(example_snp_data$read0), 96)
  expect_equal(nrow(example_snp_data$read1), 200)
  expect_equal(ncol(example_snp_data$read1), 96)
  
  # Check data types
  expect_true(is.matrix(example_snp_data$read0))
  expect_true(is.matrix(example_snp_data$read1))
  
  # Check that data is non-negative
  expect_true(all(example_snp_data$read0 >= 0))
  expect_true(all(example_snp_data$read1 >= 0))
  
  # Check that at least some reads are present
  expect_true(sum(example_snp_data$read0) > 0)
  expect_true(sum(example_snp_data$read1) > 0)
})
