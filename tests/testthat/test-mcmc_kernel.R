make_count_fixture <- function() {
  list(
    read1 = matrix(c(10, 5, 15, 8), nrow = 2),
    read0 = matrix(c(90, 95, 85, 92), nrow = 2)
  )
}

make_categorical_fixture <- function() {
  matrix(c(0, 1, 0.5, 1, 0, 0.5), nrow = 2, ncol = 3)
}

#' Sparse long-format counts: s2 has no t1 reads, s3 has no t2 (NA after grid completion).
make_count_fixture_with_na <- function() {
  data.frame(
    specimen_id = c("s1", "s1", "s1", "s2", "s3", "s3"),
    target_id = c("t1", "t1", "t2", "t2", "t1", "t1"),
    target_value = c("A", "T", "G", "G", "A", "T"),
    target_count = c(5, 3, 9, 7, 2, 6)
  )
}

expect_kernel_matches_r <- function(model_obj, seed = 42L, n_iter = 5L) {
  state_r <- run_kernel_iterations(model_obj, mcmc_kernel_r(), seed = seed, n_iter = n_iter)
  state_cpp <- run_kernel_iterations(model_obj, mcmc_kernel_cpp_ad(), seed = seed, n_iter = n_iter)
  expect_kernel_states_equal(state_cpp, state_r)
  invisible(list(r = state_r, cpp = state_cpp))
}

expect_kernel_sm_valid <- function(model_obj, seed = 42L, n_iter = 5L) {
  state <- run_kernel_iterations(model_obj, mcmc_kernel_cpp_sm(), seed = seed, n_iter = n_iter)
  expect_true(all(is.finite(state$mu)))
  expect_true(state$kplus >= 1L && state$kplus <= state$ktrunc)
  expect_equal(length(state$mu), state$ktrunc)
  expect_true(all(state$mu > 0 & state$mu <= 1))

  model_obj <- ensure_kernel_obs_cache(model_obj)
  obs <- kernel_obs_args(model_obj)
  ll_r <- model_obj$loglikelihood_matrix(state$A, state$D, model_obj)
  ll_cpp <- cpp_loglik_total(
    A = state$A,
    D = state$D,
    y = model_obj$y,
    r = model_obj$r,
    model_type = model_type_id(model_obj$name),
    llik_tab = kernel_llik_tab(model_obj),
    loglik_const = obs$loglik_const,
    obs_code = obs$obs_code
  )
  expect_equal(as.numeric(ll_cpp), as.numeric(ll_r), tolerance = 1e-8)
  invisible(state)
}

expect_kernel_cpp_reproducible <- function(model_obj, kernel, seed = 42L, n_iter = 5L) {
  state1 <- run_kernel_iterations(model_obj, kernel, seed = seed, n_iter = n_iter)
  state2 <- run_kernel_iterations(model_obj, kernel, seed = seed, n_iter = n_iter)
  expect_equal(state1$A, state2$A)
  expect_equal(state1$D, state2$D)
  expect_equal(state1$mu, state2$mu)
  expect_equal(state1$kstar, state2$kstar)
  invisible(list(first = state1, second = state2))
}

expect_kernel_states_equal <- function(state_cpp, state_r) {
  expect_equal(state_cpp$A, state_r$A)
  expect_equal(state_cpp$D, state_r$D)
  expect_equal(state_cpp$mu, state_r$mu)
  expect_equal(state_cpp$kstar, state_r$kstar)
  expect_equal(state_cpp$kplus, state_r$kplus)
  expect_equal(state_cpp$ktrunc, state_r$ktrunc)
  expect_equal(as.numeric(state_cpp$loglik), as.numeric(state_r$loglik), tolerance = 1e-6)
  expect_equal(as.numeric(state_cpp$logpost), as.numeric(state_r$logpost), tolerance = 1e-6)
}

run_kernel_iterations <- function(model_obj, kernel, seed, n_iter = 5L) {
  set.seed(seed)
  state <- model_obj$initialize_state(model_obj, threshold = 0.001)
  state <- slice_init(state, model_obj)
  model_obj$kernel <- kernel
  for (iter in seq_len(n_iter)) {
    state <- slice_iter(state, model_obj)
  }
  state
}

#' Kernel that uses compiled A/D updates with R reference s/mu (for equivalence tests)
mcmc_kernel_cpp_ad <- function() {
  list(
    name = "cpp",
    update_s = slice_update_s_r,
    update_mu = slice_update_mu_r,
    update_a = kernel_update_a_cpp,
    update_d = kernel_update_d_cpp
  )
}

#' Kernel that uses compiled s/mu updates with R reference A/D (for equivalence tests)
mcmc_kernel_cpp_sm <- function() {
  list(
    name = "cpp",
    update_s = kernel_update_s_cpp,
    update_mu = kernel_update_mu_cpp,
    update_a = slice_update_a_r,
    update_d = slice_update_d_r
  )
}

test_that("fused and sequential cpp kernels are each reproducible", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)

  expect_kernel_cpp_reproducible(model_obj, mcmc_kernel_cpp())

  kernel_seq <- mcmc_kernel_cpp()
  kernel_seq$update_iter <- NULL
  expect_kernel_cpp_reproducible(model_obj, kernel_seq)
})

test_that("compiled A/D kernel matches R reference on count models", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  expect_kernel_matches_r(model_obj)
})

test_that("compiled s/mu kernel produces valid state on count models", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  expect_kernel_sm_valid(model_obj)
})

test_that("compiled s/mu kernel produces valid state on poisson and binomial count models", {
  skip_if_not(cpp_kernel_available())

  for (model_name in c("poisson", "binomial")) {
    processed <- preprocess_data(make_count_fixture(), model_name)
    model_obj <- create_model(model_name, processed)
    expect_kernel_sm_valid(model_obj)
  }
})

test_that("full compiled kernel is reproducible and self-consistent", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  kernel_seq <- mcmc_kernel_cpp()
  kernel_seq$update_iter <- NULL
  expect_kernel_cpp_reproducible(model_obj, kernel_seq)
})

test_that("compiled kernel matches R when some genotypes are missing", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture_with_na(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  expect_true(any(is.na(model_obj$y)))
  model_obj <- ensure_kernel_obs_cache(model_obj)
  expect_true(any(model_obj$kernel_obs_code == 0L))

  expect_kernel_matches_r(model_obj)
})

test_that("compiled kernel matches R on poisson and binomial count models", {
  skip_if_not(cpp_kernel_available())

  for (model_name in c("poisson", "binomial")) {
    processed <- preprocess_data(make_count_fixture(), model_name)
    model_obj <- create_model(model_name, processed)
    expect_kernel_matches_r(model_obj)
  }
})

make_count_fixture_with_y_zero <- function() {
  list(
    read1 = matrix(c(10, 5, 0, 8), nrow = 2),
    read0 = matrix(c(90, 95, 85, 92), nrow = 2)
  )
}

test_that("compiled kernel matches R on count models with y=0 cells", {
  skip_if_not(cpp_kernel_available())

  for (model_name in c("poisson", "binomial", "negative_binomial")) {
    processed <- preprocess_data(make_count_fixture_with_y_zero(), model_name)
    model_obj <- create_model(model_name, processed)
    expect_true(any(model_obj$y == 0, na.rm = TRUE))
    model_obj <- ensure_kernel_obs_cache(model_obj)
    expect_true(any(model_obj$kernel_obs_code == 1L))
    expect_kernel_matches_r(model_obj)
  }
})

test_that("cpp total loglik matches R loglikelihood_matrix", {
  skip_if_not(cpp_kernel_available())

  cases <- list(
    list(
      label = "negative_binomial",
      model = "negative_binomial",
      processed = preprocess_data(make_count_fixture(), "negative_binomial")
    ),
    list(
      label = "poisson",
      model = "poisson",
      processed = preprocess_data(make_count_fixture(), "poisson")
    ),
    list(
      label = "binomial",
      model = "binomial",
      processed = preprocess_data(make_count_fixture(), "binomial")
    ),
    list(
      label = "categorical",
      model = "categorical",
      processed = preprocess_data(make_categorical_fixture(), "categorical")
    ),
    list(
      label = "negative_binomial_na",
      model = "negative_binomial",
      processed = preprocess_data(make_count_fixture_with_na(), "negative_binomial")
    )
  )

  for (case in cases) {
    model_obj <- create_model(case$model, case$processed)
    model_obj <- ensure_kernel_obs_cache(model_obj)
    set.seed(7L)
    state <- model_obj$initialize_state(model_obj, threshold = 0.001)
    state <- slice_init(state, model_obj)
    obs <- kernel_obs_args(model_obj)
    ll_r <- model_obj$loglikelihood_matrix(state$A, state$D, model_obj)
    ll_cpp <- cpp_loglik_total(
      A = state$A,
      D = state$D,
      y = model_obj$y,
      r = model_obj$r,
      model_type = model_type_id(model_obj$name),
      llik_tab = kernel_llik_tab(model_obj),
      loglik_const = obs$loglik_const,
      obs_code = obs$obs_code
    )
    expect_equal(as.numeric(ll_cpp), as.numeric(ll_r), tolerance = 1e-8, info = case$label)
  }
})

test_that("binomial y=0 loglik at prop boundaries matches stats::dbinom", {
  y <- 0
  n <- 10
  expect_equal(binomial_loglikelihood_vector(0, y, n), 0)
  expect_equal(binomial_loglikelihood_vector(1, y, n), -Inf)
  expect_equal(
    binomial_loglikelihood_vector(0.5, y, n),
    stats::dbinom(y, n, 0.5, log = TRUE)
  )
})

test_that("compiled A/D kernel matches R reference on categorical model", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_categorical_fixture(), "categorical")
  model_obj <- create_model("categorical", processed)

  state_r <- run_kernel_iterations(model_obj, mcmc_kernel_r(), seed = 42L, n_iter = 1L)
  state_cpp <- run_kernel_iterations(model_obj, mcmc_kernel_cpp_ad(), seed = 42L, n_iter = 1L)

  expect_equal(state_cpp$A, state_r$A)
  expect_equal(state_cpp$D, state_r$D)
  expect_equal(state_cpp$mu, state_r$mu)
  expect_equal(state_cpp$kstar, state_r$kstar)
  expect_equal(state_cpp$loglik, state_r$loglik, tolerance = 1e-6)
  expect_equal(state_cpp$logpost, state_r$logpost, tolerance = 1e-6)
})

test_that("categorical lookup table matches legacy cut/table path", {
  propvec <- c(0, 0.25, 0.5, 1, 0.99, 1.05)
  yvec <- c(0, 0.5, 1, 0, 0.5, 1)
  llik_tab <- build_categorical_llik_tab(0.05, 0.05)

  legacy <- 0
  for (idx in seq_along(propvec)) {
    y_val <- yvec[idx]
    if (is.na(y_val)) {
      next
    }
    prop_cat <- cut(propvec[idx], c(-Inf, 0, 0.99, 1.1), labels = c(0, 0.5, 1))
    y_cat <- factor(yvec[idx], levels = c(0, 0.5, 1))
    legacy <- legacy + llik_tab[as.integer(prop_cat), as.integer(y_cat)]
  }

  expect_equal(
    as.numeric(categorical_loglikelihood_vector(propvec, yvec, llik_tab = llik_tab)),
    as.numeric(legacy)
  )
})

test_that("auto kernel resolves to compiled adapter for all supported models", {
  skip_if_not(cpp_kernel_available())

  count_obj <- create_model("poisson", preprocess_data(make_count_fixture(), "poisson"))
  cat_obj <- create_model("categorical", preprocess_data(make_categorical_fixture(), "categorical"))

  expect_equal(resolve_mcmc_kernel(count_obj, backend = "auto")$name, "cpp")
  expect_equal(resolve_mcmc_kernel(cat_obj, backend = "auto")$name, "cpp")
})

test_that("kernel obs cache is built lazily and matches on-the-fly prepare", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  expect_false(kernel_obs_cache_ready(model_obj))

  model_obj <- ensure_kernel_obs_cache(model_obj)
  expect_true(kernel_obs_cache_ready(model_obj))

  on_the_fly <- cpp_build_kernel_obs_cache(
    y = model_obj$y,
    r = model_obj$r,
    model_type = model_type_id(model_obj$name)
  )

  expect_equal(model_obj$kernel_loglik_const, on_the_fly$loglik_const, tolerance = 1e-12)
  expect_equal(model_obj$kernel_obs_code, on_the_fly$obs_code)
  expect_true(any(model_obj$kernel_obs_code == 2L))

  model_obj2 <- ensure_kernel_obs_cache(model_obj)
  expect_identical(model_obj2$kernel_obs_code, model_obj$kernel_obs_code)
})

test_that("cpp slice_iter runs without explicit attach when kernel carries cache", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  model_obj <- setup_mcmc_kernel(model_obj, backend = "cpp")

  set.seed(1L)
  state <- model_obj$initialize_state(model_obj, threshold = 0.001)
  state <- slice_init(state, model_obj)
  expect_error(slice_iter(state, model_obj), NA)
  expect_error(slice_iter(state, model_obj), NA)
})

test_that("resolve_mcmc_kernel bundles obs cache on cpp kernel", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  model_obj$kernel <- resolve_mcmc_kernel(model_obj, backend = "cpp")

  expect_equal(length(model_obj$kernel$obs_code), model_obj$N * model_obj$P)
  expect_equal(length(model_obj$kernel$obs_loglik_const), model_obj$N * model_obj$P)
})

test_that("snp_slice still runs with compiled kernel when available", {
  skip_if_not(cpp_kernel_available())

  old_opt <- options(snp.slicer.mcmc_kernel = "cpp")
  on.exit(options(old_opt), add = TRUE)
  result <- snp_slice(make_count_fixture(), n_sample = 10, n_burnin = 5, verbose = FALSE)
  expect_s3_class(result, "snp_slice_results")
})

#' Mixed and single infections over a handful of strains, large enough that
#' the compiled updates change A, D and mu on most iterations.
make_mixing_count_fixture <- function() {
  set.seed(101L)
  n_host <- 16L
  n_target <- 12L
  strains <- matrix(stats::rbinom(4L * n_target, 1L, 0.4), nrow = 4L)
  allocation <- matrix(0, nrow = n_host, ncol = 4L)
  for (i in seq_len(n_host)) {
    allocation[i, sample.int(4L, if (i <= 6L) 1L else 2L)] <- 1
  }
  prop <- (allocation %*% strains) / rowSums(allocation)
  read1 <- matrix(stats::rbinom(n_host * n_target, 100L, prop), nrow = n_host)
  list(read1 = read1, read0 = 100 - read1)
}

test_that("compiled A/D kernel matches R reference on mixed infections for every count model", {
  skip_if_not(cpp_kernel_available())

  # Specimens here carry several strains, so update_a sees a range of j / n
  # proportions rather than the 0 and 1 of the two-specimen fixtures.
  for (model_name in c("poisson", "binomial", "negative_binomial")) {
    processed <- preprocess_data(make_mixing_count_fixture(), model_name)
    model_obj <- create_model(model_name, processed)
    states <- expect_kernel_matches_r(model_obj, n_iter = 5L)
    expect_gt(max(rowSums(states$cpp$A)), 1)
  }
})

test_that("compiled A/D kernel matches R reference on categorical mixed infections", {
  skip_if_not(cpp_kernel_available())

  fixture <- make_mixing_count_fixture()
  y_cat <- counts_to_categorical(fixture$read0, fixture$read1)
  expect_true(any(y_cat == 0.5, na.rm = TRUE))
  model_obj <- create_model("categorical", preprocess_data(y_cat, "categorical"))

  states <- expect_kernel_matches_r(model_obj, n_iter = 5L)
  expect_gt(max(rowSums(states$cpp$A)), 1)
})

test_that("compiled kernel entry points leave their inputs unmodified", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_mixing_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  model_obj <- setup_mcmc_kernel(model_obj, backend = "cpp")

  updates <- list(
    iter = kernel_update_iter_cpp,
    s = kernel_update_s_cpp,
    a = kernel_update_a_cpp,
    d = kernel_update_d_cpp,
    mu = kernel_update_mu_cpp
  )

  for (label in names(updates)) {
    set.seed(5L)
    state <- model_obj$initialize_state(model_obj, threshold = 0.001)
    state <- slice_init(state, model_obj)
    storage.mode(state$A) <- "double"
    storage.mode(state$D) <- "double"

    changed <- FALSE
    for (iter in seq_len(10L)) {
      # Arithmetic forces fresh copies that share no memory with the state.
      before <- list(A = state$A + 0, D = state$D + 0, mu = state$mu + 0)
      updated <- updates[[label]](state, model_obj)

      expect_identical(state$A, before$A, info = label)
      expect_identical(state$D, before$D, info = label)
      expect_identical(state$mu, before$mu, info = label)

      changed <- changed ||
        !identical(updated$A, before$A) ||
        !identical(updated$D, before$D) ||
        !identical(updated$mu, before$mu)
      # Advance with a full iteration so every update sees a moving state.
      state <- kernel_update_iter_cpp(updated, model_obj)
    }
    # Guard against a vacuous pass: the update must have changed something.
    expect_true(changed, info = label)
  }
})

test_that("compiled multinomial updates leave their inputs unmodified", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_mixing_count_fixture(), "multinomial")
  model_obj <- create_model("multinomial", processed)
  model_obj <- setup_mcmc_kernel(model_obj, backend = "cpp")
  kernel <- model_obj$kernel
  expect_equal(kernel$name, "cpp_multinomial")

  for (label in intersect(c("update_iter", "update_a", "update_d"), names(kernel))) {
    set.seed(5L)
    state <- model_obj$initialize_state(model_obj, threshold = 0.001)
    state <- slice_init(state, model_obj)
    storage.mode(state$A) <- "double"

    changed <- FALSE
    for (iter in seq_len(10L)) {
      before <- list(A = state$A + 0, D = state$D + 0L, mu = state$mu + 0)
      updated <- kernel[[label]](state, model_obj)

      expect_identical(state$A, before$A, info = label)
      expect_equal(state$D, before$D, info = label)
      expect_identical(state$mu, before$mu, info = label)

      changed <- changed ||
        !identical(updated$A, before$A) ||
        !isTRUE(all.equal(updated$D, before$D))
      state <- slice_iter(updated, model_obj)
    }
    expect_true(changed, info = label)
  }
})

test_that("stored samples and the MAP state are snapshots under the compiled kernel", {
  skip_if_not(cpp_kernel_available())

  old_opt <- options(snp.slicer.mcmc_kernel = "cpp")
  on.exit(options(old_opt), add = TRUE)

  n_burnin <- 5L
  result <- snp_slice(make_mixing_count_fixture(), model = "negative_binomial",
                      n_sample = 40, n_burnin = n_burnin, seed = 9, verbose = FALSE,
                      store_mcmc = TRUE)
  chain <- result$chains[[1]]
  samples <- chain$mcmc_samples
  model_obj <- create_model("negative_binomial", result$model_info$processed_data)

  # Each stored state must reproduce the log posterior recorded alongside it.
  recomputed <- vapply(samples, function(s) {
    model_obj$loglikelihood_matrix(s$A, s$D, model_obj) +
      logprior_a(s$A, s$mu, model_obj$alpha) +
      logprior_mu(s$mu, model_obj$alpha, model_obj$N) +
      logprior_d(s$D, model_obj$rho)
  }, numeric(1))
  stored <- vapply(samples, function(s) s$logpost, numeric(1))
  expect_equal(recomputed, stored, tolerance = 1e-8)
  expect_gt(length(unique(stored)), 1L)

  # The MAP estimate is the state at the MAP iteration, not a later one.
  map_sample <- samples[[chain$diagnostics$map_iteration - n_burnin]]
  expect_equal(map_sample$logpost, chain$diagnostics$map_logpost)
  map_matrices <- active_state_matrices(map_sample)
  expect_equal(chain$map_allocation_matrix, map_matrices$allocation_matrix)
  expect_equal(chain$map_dictionary_matrix, map_matrices$dictionary_matrix)
})

#' Matrix-list counts where one specimen has no reads at one target.
make_count_fixture_with_zero_coverage <- function() {
  fixture <- make_mixing_count_fixture()
  fixture$read1[3, 2] <- 0
  fixture$read0[3, 2] <- 0
  fixture
}

test_that("compiled kernel matches R on count models with zero-coverage cells", {
  skip_if_not(cpp_kernel_available())

  for (model_name in c("poisson", "binomial", "negative_binomial")) {
    processed <- preprocess_data(make_count_fixture_with_zero_coverage(), model_name)
    model_obj <- create_model(model_name, processed)
    expect_true(any(model_obj$r == 0))

    model_obj <- ensure_kernel_obs_cache(model_obj)
    expect_true(all(is.finite(model_obj$kernel_loglik_const)), info = model_name)

    states <- expect_kernel_matches_r(model_obj, n_iter = 3L)
    expect_true(is.finite(states$cpp$logpost), info = model_name)
  }
})

test_that("snp_slice runs with the compiled kernel on zero-coverage cells", {
  skip_if_not(cpp_kernel_available())

  old_opt <- options(snp.slicer.mcmc_kernel = "cpp")
  on.exit(options(old_opt), add = TRUE)
  result <- snp_slice(make_count_fixture_with_zero_coverage(), model = "negative_binomial",
                      n_sample = 10, n_burnin = 5, seed = 2, verbose = FALSE)
  expect_true(is.finite(result$chains[[1]]$diagnostics$final_logpost))
})

test_that("logprior_a equals the elementwise Bernoulli log prior", {
  set.seed(21L)
  A <- matrix(as.numeric(stats::runif(15 * 6) < 0.4), nrow = 15, ncol = 6)
  A[, 6] <- 0
  mu <- sort(stats::runif(6, 0.01, 0.9), decreasing = TRUE)

  elementwise <- sum(t(A) * log(mu)) + sum(t(1 - A) * log(1 - mu))
  expect_equal(logprior_a(A, mu, alpha = 2.6), elementwise, tolerance = 1e-12)
})

test_that("R reference update_d leaves fixed strains alone when kmin exceeds kstar", {
  # Every specimen is a single infection, so kmin = kstar + 1 and there is no
  # strain for update_d to resample.
  strains <- rbind(c(1, 0, 0, 1), c(0, 1, 1, 0))
  alt <- strains[rep(1:2, each = 3), ]
  processed <- preprocess_data(list(read1 = alt * 100, read0 = (1 - alt) * 100),
                               "negative_binomial")
  model_obj <- create_model("negative_binomial", processed, rho = 0.3)

  set.seed(3L)
  state <- model_obj$initialize_state(model_obj, threshold = 0.001)
  state <- slice_init(state, model_obj)
  expect_gt(state$kmin, state$kstar)

  # A dictionary entry the data contradicts: resampling strain kstar would
  # flip it back, so an unchanged D shows the strain was not touched.
  hosts <- which(state$A[, state$kstar] > 0)
  target <- which(processed$y[hosts[1], ] > 0)[1]
  state$D[state$kstar, target] <- 0

  updated <- slice_update_d_r(state, model_obj)
  expect_identical(updated$D, state$D)
})
