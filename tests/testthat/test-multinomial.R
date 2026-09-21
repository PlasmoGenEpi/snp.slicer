# Tests for the multinomial observation model: multi-allelic loading, the
# expanded dictionary layout, exact nesting of the binomial model, the Gibbs
# dictionary update, and end-to-end runs on data with more than two alleles.

# t1 triallelic (A/T/G), t2 biallelic (C/T), t3 monomorphic (G); s4 has no t2 reads.
make_multiallelic_df <- function() {
  data.frame(
    specimen_id  = c("s1", "s1", "s1", "s2", "s2", "s2", "s2", "s3", "s3", "s3", "s4", "s4"),
    target_id    = c("t1", "t2", "t3", "t1", "t1", "t2", "t3", "t1", "t2", "t3", "t1", "t3"),
    target_value = c("A",  "C",  "G",  "A",  "T",  "C",  "G",  "G",  "T",  "G",  "T",  "G"),
    target_count = c(50,   40,   30,   20,   25,   10,   12,   33,   31,   9,    15,   8)
  )
}

# Biallelic example data as a read0/read1 list, small enough for fast tests.
make_example_list <- function(n = 60, p = 12) {
  utils::data("example_snp_data", package = "snp.slicer", envir = environment())
  list(
    read0 = example_snp_data$read0[seq_len(n), seq_len(p)],
    read1 = example_snp_data$read1[seq_len(n), seq_len(p)]
  )
}

test_that("load_dataframe_multiallelic keeps every allele and orders slots by total count", {
  out <- snp.slicer:::load_dataframe_multiallelic(make_multiallelic_df())

  expect_equal(out$data_type, "read_counts_multiallelic")
  expect_equal(out$N, 4)
  expect_equal(out$P, 3)
  expect_equal(out$target_ids, c("t1", "t2", "t3"))
  expect_equal(out$specimen_ids, c("s1", "s2", "s3", "s4"))
  expect_equal(out$n_alleles, c(3L, 2L, 1L))

  # Slot 1 is the population major allele: t1 totals A 70, T 40, G 33.
  expect_equal(out$allele_labels[[1]], c("A", "T", "G"))
  expect_equal(out$allele_labels[[2]], c("C", "T"))
  expect_equal(out$allele_labels[[3]], "G")

  # Expanded layout bookkeeping.
  expect_equal(out$col_offset, c(0L, 3L, 5L))
  expect_equal(out$col_locus, c(1L, 1L, 1L, 2L, 2L, 3L))
  expect_equal(out$col_allele, c(0L, 1L, 2L, 0L, 1L, 0L))
  expect_equal(dim(out$counts_exp), c(4, 6))

  # Counts land in the right slots; missing genotype is NA across the block.
  expect_equal(unname(out$counts_exp["s2", 1:3]), c(20, 25, 0))
  expect_equal(unname(out$counts_exp["s3", 1:3]), c(0, 0, 33))
  expect_true(all(is.na(out$counts_exp["s4", 4:5])))
  expect_true(is.na(out$r["s4", "t2"]))
  expect_equal(out$r["s2", "t1"], 45)

  # Compatibility fields: y holds slot-1 counts, r0/r1 the first two labels.
  expect_equal(out$y["s2", "t1"], 20)
  expect_equal(out$r0_values, c("A", "C", "G"))
  expect_equal(out$r1_values, c("T", "T", "G"))
})

test_that("validate_input_data accepts a dataset whose only variation is triallelic for the multinomial model", {
  df <- data.frame(
    specimen_id  = c("s1", "s2", "s3", "s1"),
    target_id    = c("t1", "t1", "t1", "t2"),
    target_value = c("A",  "C",  "G",  "T"),
    target_count = c(5,    4,    3,    6)
  )
  expect_true(snp.slicer:::validate_input_data(df, "multinomial"))
  expect_error(snp.slicer:::validate_input_data(df, "negative_binomial"), "biallelic")
})

test_that("the multinomial loader does not warn about or drop polyallelic targets", {
  df <- make_multiallelic_df()
  expect_no_warning(out <- snp.slicer:::preprocess_data(df, model = "multinomial"))
  expect_true("t1" %in% out$target_ids)
  expect_warning(bi <- snp.slicer:::preprocess_data(df, model = "binomial"), "more than two alleles")
  expect_false("t1" %in% bi$target_ids)
})

test_that("expand_dictionary places one 1 per target block at the coded allele", {
  pd <- snp.slicer:::load_dataframe_multiallelic(make_multiallelic_df())
  m <- snp.slicer:::create_model("multinomial", pd, alpha = 2.6)
  D <- rbind(c(0, 0, 0), c(2, 1, 0))
  Dx <- snp.slicer:::expand_dictionary(D, m)
  expect_equal(dim(Dx), c(2, 6))
  expect_equal(Dx[1, ], c(1, 0, 0, 1, 0, 1))
  expect_equal(Dx[2, ], c(0, 0, 1, 0, 1, 1))
  expect_equal(dim(snp.slicer:::expand_dictionary(D[0, , drop = FALSE], m)), c(0, 6))
})

test_that("dict_prior options resolve to per-target probability vectors", {
  pd <- snp.slicer:::load_dataframe_multiallelic(make_multiallelic_df())
  emp <- snp.slicer:::resolve_dict_prior("empirical", pd)
  expect_equal(lengths(emp), c(3L, 2L, 1L))
  expect_equal(unname(emp[[1]]), (c(70, 40, 33) + 1) / sum(c(70, 40, 33) + 1))
  uni <- snp.slicer:::resolve_dict_prior("uniform", pd)
  expect_equal(unname(uni[[1]]), rep(1 / 3, 3))
  custom <- snp.slicer:::resolve_dict_prior(list(c(2, 1, 1), c(1, 1), 5), pd)
  expect_equal(custom[[1]], c(0.5, 0.25, 0.25))
  expect_equal(custom[[3]], 1)
  expect_error(snp.slicer:::resolve_dict_prior(list(c(1, 1), c(1, 1), 1), pd), "3 non-negative weights")
  expect_error(snp.slicer:::resolve_dict_prior("flat", pd), "dict_prior must be")
})

test_that("multinomial likelihood equals the binomial likelihood on two-allele data", {
  dat <- make_example_list()
  mb <- snp.slicer:::create_model("binomial", snp.slicer:::preprocess_data(dat, "binomial"), alpha = 2.6)
  mm <- snp.slicer:::create_model("multinomial", snp.slicer:::preprocess_data(dat, "multinomial"), alpha = 2.6)

  set.seed(7)
  sb <- mb$initialize_state(mb, 0.001)
  # Binomial: y = read1, D == 1 carries the read1 allele. List input keeps read0
  # in slot 1, so multinomial code 1 is the read1 allele too: same D.
  expect_equal(
    mm$loglikelihood_matrix(sb$A, sb$D, mm),
    mb$loglikelihood_matrix(sb$A, sb$D, mb),
    tolerance = 1e-10
  )

  # The kernel (vector) form drops the multinomial coefficient, which for two
  # alleles is exactly the binomial coefficient dbinom() includes.
  i <- sb$mixed[1]
  prop <- as.vector((sb$A[i, , drop = FALSE] %*% sb$D) / sum(sb$A[i, ]))
  propx <- as.vector((sb$A[i, , drop = FALSE] %*% snp.slicer:::expand_dictionary(sb$D, mm)) / sum(sb$A[i, ]))
  expect_equal(
    mb$loglikelihood_vector(prop, mb$y[i, ], mb$r[i, ]) - sum(lchoose(mb$r[i, ], mb$y[i, ]), na.rm = TRUE),
    mm$loglikelihood_vector(propx, mm$counts_exp[i, ]),
    tolerance = 1e-10
  )
})

test_that("Gibbs cell log weights match brute-force full posterior differences", {
  dat <- make_example_list()
  mm <- snp.slicer:::create_model("multinomial", snp.slicer:::preprocess_data(dat, "multinomial"), alpha = 2.6)
  set.seed(7)
  st <- mm$initialize_state(mm, 0.001)
  A <- st$A
  D <- st$D
  active <- which(colSums(A) > 0)
  k <- active[active >= st$kmin][1]
  p <- 3L
  hosts <- which(A[, k] > 0)
  cols <- mm$col_offset[p] + seq_len(mm$n_alleles[p])

  M <- A %*% snp.slicer:::expand_dictionary(D, mm)
  base <- M[hosts, cols, drop = FALSE]
  base[, D[k, p] + 1] <- base[, D[k, p] + 1] - A[hosts, k]
  lw <- snp.slicer:::multinomial_cell_logweights(
    base, mm$counts_exp[hosts, cols, drop = FALSE], A[hosts, k],
    mm$log_dict_prior_pad[seq_along(cols), p]
  )
  brute <- vapply(seq_along(cols) - 1, function(m) {
    D2 <- D
    D2[k, p] <- m
    mm$loglikelihood_matrix(A, D2, mm) + snp.slicer:::multinomial_logprior_d(D2, mm)
  }, numeric(1))

  # Equal up to the constant shared by every allele option.
  expect_equal(diff(lw - brute), 0, tolerance = 1e-8)
})

test_that("multinomial_logprior_d sums the per-target log prior of each code", {
  pd <- snp.slicer:::load_dataframe_multiallelic(make_multiallelic_df())
  m <- snp.slicer:::create_model("multinomial", pd, alpha = 2.6, dict_prior = "uniform")
  D <- rbind(c(0, 1, 0), c(2, 0, 0))
  expect_equal(snp.slicer:::multinomial_logprior_d(D, m), 2 * (log(1 / 3) + log(1 / 2) + 0))
  expect_equal(snp.slicer:::multinomial_logprior_d(D[0, , drop = FALSE], m), 0)
})

test_that("multinomial initialization is valid and repairs unexplained alleles", {
  pd <- snp.slicer:::load_dataframe_multiallelic(make_multiallelic_df())
  m <- snp.slicer:::create_model("multinomial", pd, alpha = 2.6)
  set.seed(1)
  st <- m$initialize_state(m, 0.001)

  expect_true(is.finite(st$loglik))
  # Every strain has exactly one allele code per target, within range.
  for (p in seq_len(m$P)) {
    expect_true(all(st$D[, p] %in% (seq_len(m$n_alleles[p]) - 1)))
  }
  # s2 shows A and T at t1: it must carry strains giving both alleles.
  Dx <- snp.slicer:::expand_dictionary(st$D, m)
  carried <- (st$A %*% Dx) > 0
  expect_true(carried[2, 1] && carried[2, 2])
  # Single infections (s1, s3, s4) carry exactly one strain; s2 is the only mixed one.
  expect_equal(st$mixed, 2L)
  expect_equal(unname(rowSums(st$A)[c(1, 3, 4)]), c(1, 1, 1))
})

test_that("snp_slice runs end-to-end with the multinomial model on triallelic data", {
  df <- make_multiallelic_df()
  res <- snp_slice(df, model = "multinomial", n_sample = 30, n_burnin = 30, seed = 3, verbose = FALSE)

  expect_s3_class(res, "snp_slice_results")
  expect_equal(res$model_info$model, "multinomial")
  expect_equal(res$model_info$data_type, "read_counts_multiallelic")
  chain <- get_chain(res)
  expect_true(is.finite(chain$diagnostics$map_logpost))
  expect_true(all(chain$map_dictionary_matrix[, 1] %in% 0:2))
  expect_true(all(chain$map_dictionary_matrix[, 2] %in% 0:1))
  expect_true(all(chain$map_dictionary_matrix[, 3] == 0))

  af <- calculate_allele_frequencies(res, c("t1", "t2", "t3"), estimate = "map")
  expect_equal(sum(af$frequency), 1)
  # The three single infections fix A|C|G, G|T|G and T|C|G; s2 (A/T at t1, C at t2)
  # must be explained by A|C|G plus T|C|G. So those haplotypes carry all the mass.
  present <- af$allele[af$frequency > 0]
  expect_setequal(present, c("A|C|G", "T|C|G", "G|T|G"))
  expect_equal(af$count[af$allele == "G|T|G"], 1)
})

test_that("snp_slice_multinomial accepts read0/read1 lists and the posterior estimator", {
  dat <- make_example_list(n = 30, p = 8)
  res <- snp_slice_multinomial(dat, n_sample = 20, n_burnin = 10, seed = 5,
                               verbose = FALSE, store_mcmc = TRUE)
  expect_equal(res$model_info$model, "multinomial")
  expect_equal(res$model_info$processed_data$allele_labels[[1]], c("ref", "alt"))
  af <- calculate_allele_frequencies(res, 1:3, estimate = "posterior", n_samples = 10)
  expect_equal(sum(af$frequency), 1)
  expect_true(all(grepl("^(ref|alt)\\|(ref|alt)\\|(ref|alt)$", af$allele)))
})

test_that("the multinomial model resolves to its fused compiled kernel", {
  skip_if_not(snp.slicer:::cpp_kernel_available(), "compiled kernel not available")
  pd <- snp.slicer:::load_dataframe_multiallelic(make_multiallelic_df())
  m <- snp.slicer:::create_model("multinomial", pd, alpha = 2.6)
  for (backend in c("cpp", "auto")) {
    k <- snp.slicer:::resolve_mcmc_kernel(m, backend = backend)
    expect_equal(k$name, "cpp_multinomial")
    expect_identical(k$update_iter, snp.slicer:::multinomial_update_iter_cpp)
    expect_identical(k$update_a, snp.slicer:::multinomial_update_a_cpp)
    expect_identical(k$update_d, snp.slicer:::multinomial_update_d_cpp)
    expect_identical(k$loglik, snp.slicer:::multinomial_loglik_cpp)
  }
  kr <- snp.slicer:::resolve_mcmc_kernel(m, backend = "r")
  expect_equal(kr$name, "r")
  expect_identical(kr$update_a, snp.slicer:::multinomial_update_a_r)
  expect_identical(kr$update_d, snp.slicer:::multinomial_update_d_r)
  expect_null(kr$loglik)
})

# The compiled allocation update consumes one uniform draw per (specimen,
# strain) cell exactly as the R reference does, so the two kernels must walk
# the same path from the same seed.
run_multinomial_kernel <- function(model_obj, kernel, seed, n_iter) {
  set.seed(seed)
  state <- model_obj$initialize_state(model_obj, threshold = 0.001)
  state <- snp.slicer:::slice_init(state, model_obj)
  model_obj$kernel <- kernel
  for (iter in seq_len(n_iter)) {
    state <- snp.slicer:::slice_iter(state, model_obj)
  }
  state
}

multinomial_fixtures <- function() {
  list(
    triallelic = snp.slicer:::preprocess_data(make_multiallelic_df(), "multinomial"),
    example = snp.slicer:::preprocess_data(make_example_list(n = 60, p = 12), "multinomial")
  )
}

test_that("compiled multinomial A and D updates match the R reference from the same seed", {
  skip_if_not(snp.slicer:::cpp_kernel_available(), "compiled kernel not available")
  # fused = FALSE keeps the R slice-variable and stick-breaking updates, whose
  # grid sampler draws differently in C++, so the comparison isolates the
  # allocation and dictionary updates.
  fixtures <- multinomial_fixtures()
  for (nm in names(fixtures)) {
    m <- snp.slicer:::create_model("multinomial", fixtures[[nm]], alpha = 2.6)
    for (seed in c(42L, 7L)) {
      r_state <- run_multinomial_kernel(m, snp.slicer:::mcmc_kernel_r(m), seed = seed, n_iter = 10L)
      c_state <- run_multinomial_kernel(m, snp.slicer:::mcmc_kernel_multinomial_cpp(m, fused = FALSE),
                                        seed = seed, n_iter = 10L)
      expect_identical(c_state$A, r_state$A, info = nm)
      expect_identical(c_state$D, r_state$D, info = nm)
      expect_equal(c_state$kstar, r_state$kstar, info = nm)
      expect_equal(c_state$loglik, r_state$loglik, info = nm)
      expect_equal(c_state$logpost, r_state$logpost, info = nm)
    }
  }
})

test_that("compiled multinomial likelihood matches the R routine on valid and impossible states", {
  skip_if_not(snp.slicer:::cpp_kernel_available(), "compiled kernel not available")
  m <- snp.slicer:::create_model("multinomial", multinomial_fixtures()$example, alpha = 2.6)
  set.seed(9)
  st <- m$initialize_state(m, 0.001)
  expect_equal(snp.slicer:::multinomial_loglik_cpp(st, m), m$loglikelihood_matrix(st$A, st$D, m), tolerance = 1e-10)
  # Remove every strain but one from a mixed specimen: some allele it shows loses its carrier.
  i <- st$mixed[1]
  A_bad <- st$A
  keep <- which(A_bad[i, ] > 0)[1]
  A_bad[i, ] <- 0
  A_bad[i, keep] <- 1
  bad <- list(A = A_bad, D = st$D)
  expect_equal(snp.slicer:::multinomial_loglik_cpp(bad, m), m$loglikelihood_matrix(A_bad, st$D, m))
})

test_that("the fused compiled multinomial iteration is reproducible and self-consistent", {
  skip_if_not(snp.slicer:::cpp_kernel_available(), "compiled kernel not available")
  fixtures <- multinomial_fixtures()
  for (nm in names(fixtures)) {
    m <- snp.slicer:::create_model("multinomial", fixtures[[nm]], alpha = 2.6)
    k <- snp.slicer:::mcmc_kernel_multinomial_cpp(m)
    s1 <- run_multinomial_kernel(m, k, seed = 42L, n_iter = 10L)
    s2 <- run_multinomial_kernel(m, k, seed = 42L, n_iter = 10L)
    expect_identical(s1$A, s2$A, info = nm)
    expect_identical(s1$D, s2$D, info = nm)
    expect_identical(s1$mu, s2$mu, info = nm)
    # Shapes and bookkeeping stay coherent as the stick-breaking grows the state.
    expect_equal(ncol(s1$A), nrow(s1$D), info = nm)
    expect_equal(length(s1$mu), s1$ktrunc, info = nm)
    expect_true(s1$kplus >= 1L && s1$kplus <= s1$ktrunc, info = nm)
    expect_true(all(s1$mu > 0 & s1$mu <= 1), info = nm)
    for (p in seq_len(m$P)) {
      expect_true(all(s1$D[, p] %in% (seq_len(m$n_alleles[p]) - 1)), info = nm)
    }
    expect_true(all(rowSums(s1$A) >= 1), info = nm)
    # The compiled likelihood the adapter reports equals the R routine on the final state.
    expect_equal(s1$loglik, m$loglikelihood_matrix(s1$A, s1$D, m), tolerance = 1e-10, info = nm)
    expect_true(is.finite(s1$logpost), info = nm)
  }
})

test_that("the compiled multinomial cell likelihood matches the R kernel term, including impossible states", {
  skip_if_not(snp.slicer:::cpp_kernel_available(), "compiled kernel not available")
  pd <- snp.slicer:::preprocess_data(make_multiallelic_df(), "multinomial")
  m <- snp.slicer:::create_model("multinomial", pd, alpha = 2.6)
  set.seed(3)
  st <- m$initialize_state(m, 0.001)
  Dx <- snp.slicer:::expand_dictionary(st$D, m)
  cache <- snp.slicer:::cpp_build_kernel_obs_cache(m$counts_exp, m$counts_exp, snp.slicer:::model_type_id("multinomial"))
  # Valid state: compiled total equals the R kernel term (no coefficients on either side).
  ll_cpp <- snp.slicer:::cpp_loglik_total(st$A, Dx, m$counts_exp, m$counts_exp, snp.slicer:::model_type_id("multinomial"),
                                          snp.slicer:::NULL_LLIK_TAB, cache$loglik_const, cache$obs_code)
  ll_r <- m$loglikelihood_matrix(st$A, st$D, m) - m$loglik_const
  expect_equal(as.numeric(ll_cpp), ll_r, tolerance = 1e-10)
  # Impossible state: drop the only strain carrying an allele a specimen shows reads for.
  A_bad <- st$A
  A_bad[2, ] <- 0
  A_bad[2, which(st$D[, 1] == 0)[1]] <- 1   # s2 shows A and T at t1; keep only an A carrier
  ll_bad <- snp.slicer:::cpp_loglik_total(A_bad, Dx, m$counts_exp, m$counts_exp, snp.slicer:::model_type_id("multinomial"),
                                          snp.slicer:::NULL_LLIK_TAB, cache$loglik_const, cache$obs_code)
  expect_equal(as.numeric(ll_bad), -Inf)
  expect_equal(m$loglikelihood_matrix(A_bad, st$D, m), -Inf)
})

test_that("multinomial dictionary codes are labelled with the per-target allele labels", {
  A <- matrix(c(1, 1, 0, 1), nrow = 2)               # s1: strain 1; s2: strains 1 and 2
  D <- rbind(c(2, 0), c(1, 1))                        # codes into the label lists below
  result <- structure(list(
    map_allocation_matrix = A, map_dictionary_matrix = D,
    final_allocation_matrix = A, final_dictionary_matrix = D,
    model_info = list(
      model = "multinomial", N = 2, P = 2, data_type = "read_counts_multiallelic",
      processed_data = list(
        target_ids = c("t1", "t2"),
        r0_values = c("A", "C"), r1_values = c("T", "T"),
        allele_labels = list(c("A", "T", "G"), c("C", "T"))
      )
    )
  ), class = "snp_slice_results")
  af <- calculate_allele_frequencies(result, c(1L, 2L), estimate = "map")
  expect_equal(af$count[af$allele == "G|C"], 2)
  expect_equal(af$count[af$allele == "T|T"], 1)
  expect_equal(sum(af$count), 3)
  expect_equal(nrow(af), 6)
})
