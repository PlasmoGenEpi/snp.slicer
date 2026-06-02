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
  expect_equal(state_cpp$A, state_r$A)
  expect_equal(state_cpp$D, state_r$D)
  expect_equal(state_cpp$mu, state_r$mu)
  expect_equal(state_cpp$kstar, state_r$kstar)
  expect_equal(state_cpp$kplus, state_r$kplus)
  expect_equal(state_cpp$ktrunc, state_r$ktrunc)
  expect_equal(as.numeric(state_cpp$loglik), as.numeric(state_r$loglik), tolerance = 1e-6)
  expect_equal(as.numeric(state_cpp$logpost), as.numeric(state_r$logpost), tolerance = 1e-6)
  invisible(list(r = state_r, cpp = state_cpp))
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

test_that("fused cpp_slice_iter matches sequential cpp kernel updates", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)

  run_iters <- function(kernel, seed = 42L) {
    set.seed(seed)
    state <- model_obj$initialize_state(model_obj, threshold = 0.001)
    state <- slice_init(state, model_obj)
    model_obj$kernel <- kernel
    for (iter in seq_len(5L)) {
      state <- slice_iter(state, model_obj)
    }
    state
  }

  kernel_seq <- mcmc_kernel_cpp()
  kernel_seq$update_iter <- NULL

  state_fused <- run_iters(mcmc_kernel_cpp())
  state_seq <- run_iters(kernel_seq)

  expect_equal(state_fused$A, state_seq$A)
  expect_equal(state_fused$D, state_seq$D)
  expect_equal(state_fused$mu, state_seq$mu)
  expect_equal(state_fused$kstar, state_seq$kstar)
  expect_equal(state_fused$kplus, state_seq$kplus)
  expect_equal(state_fused$ktrunc, state_seq$ktrunc)
})

test_that("compiled A/D kernel matches R reference on count models", {
  skip_if_not(cpp_kernel_available())

  processed <- preprocess_data(make_count_fixture(), "negative_binomial")
  model_obj <- create_model("negative_binomial", processed)
  expect_kernel_matches_r(model_obj)
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
  result <- snp_slice(make_count_fixture(), n_mcmc = 25, verbose = FALSE)
  expect_s3_class(result, "snp_slice_results")
})
