# Tests for the burn-in contract: the chain runs n_burnin + n_sample iterations
# and only the post-burn-in iterations are retained, so every downstream
# consumer sees a chain that is post-burn-in by construction.

make_burnin_data <- function(n_hosts = 6, n_snps = 5, seed = 42) {
  set.seed(seed)
  read1 <- matrix(stats::rpois(n_hosts * n_snps, 20), nrow = n_hosts, ncol = n_snps)
  read0 <- matrix(stats::rpois(n_hosts * n_snps, 20), nrow = n_hosts, ncol = n_snps)
  list(read1 = read1, read0 = read0)
}

test_that("burn-in iterations are run but not retained", {
  result <- snp_slice(make_burnin_data(), n_sample = 20, n_burnin = 10,
                      store_mcmc = TRUE, verbose = FALSE)

  expect_equal(length(get_chain(result)$mcmc_samples), 20)
  expect_equal(get_chain(result)$convergence$iterations_run, 30)
  expect_equal(get_chain(result)$convergence$samples_retained, 20)
  expect_equal(result$parameters$n_sample, 20)
  expect_equal(result$parameters$n_burnin, 10)
})

test_that("retained samples carry their true iteration number and ktrunc", {
  result <- snp_slice(make_burnin_data(), n_sample = 15, n_burnin = 8,
                      store_mcmc = TRUE, verbose = FALSE)

  iterations <- vapply(get_chain(result)$mcmc_samples, function(s) s$iteration, numeric(1))
  expect_equal(iterations, 9:23)
  expect_true(all(vapply(get_chain(result)$mcmc_samples, function(s) !is.null(s$ktrunc), logical(1))))
})

test_that("MAP state is selected from post-burn-in iterations only", {
  result <- snp_slice(make_burnin_data(), n_sample = 20, n_burnin = 10,
                      store_mcmc = TRUE, verbose = FALSE)

  expect_gt(get_chain(result)$diagnostics$map_iteration, 10)
})

test_that("MCMC consumers agree on the sample pool they use", {
  result <- snp_slice(make_burnin_data(), n_sample = 20, n_burnin = 10,
                      store_mcmc = TRUE, verbose = FALSE)

  coi <- calculate_individual_coi(result, estimate = "posterior", n_samples = 20)
  expect_equal(nrow(coi), 6)
  expect_false(any(is.na(coi$coi_estimate)))

  freqs <- calculate_allele_frequencies(result, c(1, 2), estimate = "posterior", n_samples = 20)
  expect_equal(unique(freqs$n_samples), 20)

  draws <- as_draws_snp_slice(result, "logpost")
  expect_equal(posterior::niterations(draws), 20)
})

test_that("additional_burnin trims the retained chain further", {
  result <- snp_slice(make_burnin_data(), n_sample = 20, n_burnin = 10,
                      store_mcmc = TRUE, verbose = FALSE)

  draws <- as_draws_snp_slice(result, "logpost", additional_burnin = 15)
  expect_equal(posterior::niterations(draws), 5)

  expect_error(
    as_draws_snp_slice(result, "logpost", additional_burnin = 20),
    "additional_burnin"
  )
})

test_that("plot_convergence uses true iteration numbers", {
  skip_if_not_installed("ggplot2")
  result <- snp_slice(make_burnin_data(), n_sample = 20, n_burnin = 10,
                      store_mcmc = TRUE, verbose = FALSE)

  p <- plot_convergence(result, "logpost")
  expect_equal(min(p$data$iteration), 11)
  expect_equal(max(p$data$iteration), 30)
})
