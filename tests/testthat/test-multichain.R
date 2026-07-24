make_chain_data <- function(n = 15, p = 10, seed = 7) {
  set.seed(seed)
  list(
    read1 = matrix(stats::rpois(n * p, 20), nrow = n, ncol = p),
    read0 = matrix(stats::rpois(n * p, 20), nrow = n, ncol = p)
  )
}

run_chains_fixture <- function(n_chains = 3, n_cores = 1, seed = 42) {
  snp_slice(make_chain_data(), n_sample = 10, n_burnin = 5,
            n_chains = n_chains, n_cores = n_cores, seed = seed,
            store_mcmc = TRUE, verbose = FALSE)
}

test_that("n_chains runs the requested number of chains with distinct seeds", {
  result <- run_chains_fixture(n_chains = 3)

  expect_length(result$chains, 3)
  comparison <- compare_chains(result)
  expect_equal(nrow(comparison), 3)
  expect_equal(comparison$chain, 1:3)
  expect_equal(anyDuplicated(comparison$seed), 0)
  expect_equal(sum(comparison$best), 1)
})

test_that("best_chain is the chain with the highest MAP logpost", {
  result <- run_chains_fixture(n_chains = 3)

  comparison <- compare_chains(result)
  best <- which.max(comparison$map_logpost)
  expect_equal(result$best_chain, best)
  expect_identical(get_chain(result)$allocation_matrix,
                   result$chains[[best]]$allocation_matrix)
})

test_that("multi-chain runs are reproducible from a single seed", {
  a <- run_chains_fixture(n_chains = 2, seed = 99)
  b <- run_chains_fixture(n_chains = 2, seed = 99)

  expect_equal(compare_chains(a), compare_chains(b))
  expect_identical(get_chain(a)$allocation_matrix, get_chain(b)$allocation_matrix)
})

test_that("parallel dispatch gives the same chains as sequential dispatch", {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2, "requires more than one core")

  serial <- run_chains_fixture(n_chains = 2, n_cores = 1, seed = 5)
  parallel_run <- run_chains_fixture(n_chains = 2, n_cores = 2, seed = 5)

  expect_equal(compare_chains(serial)$map_logpost,
               compare_chains(parallel_run)$map_logpost)
  expect_identical(get_chain(serial)$allocation_matrix, get_chain(parallel_run)$allocation_matrix)
})

test_that("get_chain returns a usable single-chain results object", {
  result <- run_chains_fixture(n_chains = 2)

  chain2 <- get_chain(result, 2)
  expect_s3_class(chain2, "snp_slice_results")
  expect_identical(chain2$allocation_matrix, result$chains[[2]]$allocation_matrix)
  expect_identical(chain2$mcmc_samples, result$chains[[2]]$mcmc_samples)
  expect_equal(chain2$model_info$model, result$model_info$model)

  # Downstream diagnostics work unchanged on a single chain
  coi <- calculate_individual_coi(chain2, use_map = FALSE, n_samples = 5)
  expect_equal(nrow(coi), nrow(chain2$allocation_matrix))
  expect_true(is.data.frame(convergence_diagnostics(chain2, pars = "logpost")))

  # Default is the best chain
  expect_identical(get_chain(result)$allocation_matrix,
                   result$chains[[result$best_chain]]$allocation_matrix)
  expect_error(get_chain(result, 5), "chain must be an integer")
})

test_that("diagnostics take a chain argument, defaulting to the best chain", {
  result <- run_chains_fixture(n_chains = 3)
  other <- setdiff(seq_along(result$chains), result$best_chain)[1]

  # NULL (default) means the top-level, i.e. best-chain, estimates
  expect_identical(extract_strains(result), extract_strains(result, chain = result$best_chain))
  expect_identical(extract_allocations(result)$allocation_matrix,
                   get_chain(result)$allocation_matrix)

  # An explicit index selects that chain instead
  expect_identical(extract_strains(result, chain = other)$dictionary,
                   result$chains[[other]]$dictionary_matrix)
  expect_identical(extract_allocations(result, chain = other)$allocation_matrix,
                   result$chains[[other]]$allocation_matrix)

  coi <- calculate_individual_coi(result, chain = other)
  expect_equal(coi$coi_estimate, unname(rowSums(result$chains[[other]]$allocation_matrix)))

  freqs <- calculate_allele_frequencies(result, c(1, 2), chain = other)
  expect_true(all(freqs$frequency >= 0))
})

test_that("verbose multi-chain runs report progress without error", {
  expect_output(
    snp_slice(make_chain_data(), n_sample = 5, n_burnin = 3,
              n_chains = 2, seed = 11, verbose = TRUE),
    "Best chain:"
  )
})

test_that("chain settings are validated", {
  data <- make_chain_data()
  expect_error(snp_slice(data, n_chains = 0), "n_chains must be a positive integer")
  expect_error(snp_slice(data, n_chains = 1.5), "n_chains must be a positive integer")
  expect_error(snp_slice(data, n_cores = 0), "n_cores must be a positive integer")
  expect_error(snp_slice(data, n_cores = -2), "n_cores must be a positive integer")
})
