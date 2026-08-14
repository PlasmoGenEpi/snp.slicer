# Tests for R-hat / ESS convergence diagnostics (posterior-backed)

make_conv_data <- function() {
  set.seed(1)
  n <- 12
  p <- 8
  read1 <- matrix(rpois(n * p, 5), nrow = n)
  read0 <- matrix(rpois(n * p, 5), nrow = n)
  list(read1 = read1, read0 = read0)
}

run_conv_fixture <- function(n_chains = 3, seed = 1) {
  snp_slice(make_conv_data(), model = "negative_binomial",
            n_sample = 30, n_burnin = 10, n_chains = n_chains,
            seed = seed, store_mcmc = TRUE, verbose = FALSE)
}

test_that("convergence_diagnostics pools all chains into one row per parameter", {
  result <- run_conv_fixture(n_chains = 3)

  diag <- convergence_diagnostics(result)
  expect_s3_class(diag, "data.frame")
  expect_setequal(diag$variable, c("logpost", "n_strains", "kstar", "ktrunc"))
  expect_equal(nrow(diag), 4)
  expect_equal(names(diag),
               c("variable", "mean", "median", "sd", "q5", "q95",
                 "rhat", "ess_bulk", "ess_tail"))
})

test_that("as_draws_snp_slice returns an iteration x chain x variable array", {
  result <- run_conv_fixture(n_chains = 3)

  draws <- as_draws_snp_slice(result, pars = c("logpost", "n_strains"))
  expect_s3_class(draws, "draws_array")
  expect_equal(posterior::nchains(draws), 3)
  expect_equal(posterior::niterations(draws), 30)
  expect_setequal(posterior::variables(draws), c("logpost", "n_strains"))
})

test_that("single-chain runs still yield a (split-)R-hat", {
  result <- run_conv_fixture(n_chains = 1)

  diag <- convergence_diagnostics(result, pars = "logpost")
  expect_equal(nrow(diag), 1)
  # posterior splits the one chain in half, so rhat is defined (finite here as
  # logpost moves)
  expect_true(is.finite(diag$rhat))
})

test_that("coi expands to one row/variable per host", {
  result <- run_conv_fixture(n_chains = 2)
  n_hosts <- nrow(get_chain(result)$map_allocation_matrix)

  diag <- convergence_diagnostics(result, pars = "coi")
  expect_equal(nrow(diag), n_hosts)
  expect_true(all(grepl("^coi\\[[0-9]+\\]$", diag$variable)))
})

test_that("ragged chains warn (naming the shortest length) and truncate", {
  result <- run_conv_fixture(n_chains = 2)

  # Simulate gap early-stopping by dropping samples from one chain
  result$chains[[2]]$mcmc_samples <- result$chains[[2]]$mcmc_samples[1:18]

  expect_warning(
    draws <- as_draws_snp_slice(result, pars = "logpost"),
    "shortest = 18"
  )
  expect_equal(posterior::niterations(draws), 18)
})

test_that("missing MCMC samples produce a clear error", {
  result <- structure(
    list(chains = list(list(mcmc_samples = NULL)), best_chain = 1),
    class = "snp_slice_results"
  )
  expect_error(convergence_diagnostics(result), "MCMC samples not stored")
})
