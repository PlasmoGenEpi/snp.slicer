# Tests for categorical model with long-format data.frame (counts -> 0/0.5/1)

test_that("counts_to_categorical maps read counts to categories", {
  # matrix() fills by column: [10,0,3,0] -> col1 (10,0), col2 (3,0)
  read0 <- matrix(c(10, 0, 3, 0), nrow = 2, ncol = 2, dimnames = list(c("s1", "s2"), c("t1", "t2")))
  read1 <- matrix(c(0, 5, 4, 0), nrow = 2, ncol = 2, dimnames = list(c("s1", "s2"), c("t1", "t2")))
  y <- snp.slicer:::counts_to_categorical(read0, read1)
  expect_equal(y["s1", "t1"], 0)    # ref only (10, 0)
  expect_equal(y["s2", "t1"], 1)    # alt only (0, 5)
  expect_equal(y["s1", "t2"], 0.5)  # both (3, 4)
  expect_true(is.na(y["s2", "t2"])) # zero total
  expect_equal(dim(y), c(2, 2))
  expect_equal(dimnames(y), dimnames(read0))
})

test_that("counts_to_categorical errors on invalid input", {
  read0 <- matrix(c(1, 2), nrow = 1)
  read1 <- matrix(c(1, 2, 3), nrow = 1)
  expect_error(snp.slicer:::counts_to_categorical(read0, read1), "identical dimensions")
  expect_error(snp.slicer:::counts_to_categorical(matrix(-1, 1, 1), matrix(0, 1, 1)), "cannot be negative")
})

test_that("load_dataframe_categorical produces categorical processed_data from long df", {
  # Long format: 2 specimens, 2 targets, ref/alt counts -> 0, 1, 0.5, NA
  df <- data.frame(
    specimen_id = c("s1", "s1", "s1", "s1", "s2", "s2", "s2", "s2"),
    target_id   = c("t1", "t1", "t2", "t2", "t1", "t1", "t2", "t2"),
    target_value = c("ref", "alt", "ref", "alt", "ref", "alt", "ref", "alt"),
    target_count = c(10, 0, 0, 5, 3, 4, 0, 0)
  )
  out <- snp.slicer:::load_dataframe_categorical(df, model = "categorical")
  expect_equal(out$data_type, "categorical")
  expect_null(out$r)
  expect_equal(out$N, 2)
  expect_equal(out$P, 2)
  expect_equal(out$y["s1", "t1"], 0)
  expect_equal(out$y["s1", "t2"], 1)
  expect_equal(out$y["s2", "t1"], 0.5)
  expect_true(is.na(out$y["s2", "t2"]))
})

test_that("preprocess_data routes categorical data.frame through load_dataframe_categorical", {
  df <- data.frame(
    specimen_id = c("s1", "s1", "s2", "s2"),
    target_id   = c("t1", "t2", "t1", "t2"),
    target_value = c("ref", "ref", "alt", "alt"),
    target_count = c(5, 0, 0, 3)
  )
  # One value per target: ref for t1/t2 on s1, alt for t1/t2 on s2; complete will fill 0 for missing
  out <- snp.slicer:::preprocess_data(df, model = "categorical")
  expect_equal(out$data_type, "categorical")
  expect_equal(out$model, "categorical")
  expect_true(is.matrix(out$y))
  expect_equal(nrow(out$y), 2)
  expect_equal(ncol(out$y), 2)
})

test_that("snp_slice_categorical runs with long-format data.frame", {
  df <- data.frame(
    specimen_id = c("s1", "s1", "s2", "s2"),
    target_id   = c("t1", "t2", "t1", "t2"),
    target_value = c("ref", "ref", "alt", "alt"),
    target_count = c(5, 0, 0, 3)
  )
  result <- snp_slice_categorical(df, n_sample = 25, verbose = FALSE)
  expect_s3_class(result, "snp_slice_results")
  expect_equal(result$model_info$model, "categorical")
  expect_equal(result$model_info$data_type, "categorical")
})

test_that("categorical data.frame with zero total count yields NA and does not break run", {
  df <- data.frame(
    specimen_id = c("s1", "s1", "s2", "s2"),
    target_id   = c("t1", "t2", "t1", "t2"),
    target_value = c("ref", "ref", "alt", "alt"),
    target_count = c(1, 0, 0, 0)  # s2,t2 has total 0 -> NA
  )
  out <- snp.slicer:::preprocess_data(df, model = "categorical")
  expect_true(is.na(out$y["s2", "t2"]))
  expect_no_error(
    snp_slice_categorical(df, n_sample = 25, verbose = FALSE)
  )
})

test_that("validate_input_data rejects categorical data.frame with missing columns", {
  df <- data.frame(specimen_id = "s1", target_id = "t1", target_value = "ref")  # no target_count
  expect_error(
    snp.slicer:::validate_input_data(df, "categorical"),
    "missing"
  )
})

test_that("validate_input_data rejects categorical data.frame with negative target_count", {
  df <- data.frame(
    specimen_id = "s1", target_id = "t1", target_value = c("ref", "alt"),
    target_count = c(5, -1)
  )
  expect_error(
    snp.slicer:::validate_input_data(df, "categorical"),
    "cannot be negative"
  )
})

# Regression: the exception detector must classify proportions with the same
# bin edges as the likelihood. llik_tab is -Inf at [prop bin 3, y == 0], and
# bin 3 is `prop > 0.99` -- not `prop == 1`. Detecting with exact equality left
# any proportion in (0.99, 1) unfixable, so initialization aborted after 10
# no-op repair passes. Seen in the wild at 124/125 strains = 0.992.
test_that("resolve_exceptions catches proportions in (0.99, 1), not just == 1", {
  n_strain <- 125
  # one specimen, one locus observed as reference (y == 0)
  y <- matrix(0, nrow = 1, ncol = 1, dimnames = list("s1", "t1"))
  # the specimen carries all 125 strains, of which 124 have the alt allele,
  # so the proportion is 124/125 = 0.992 -- inside (0.99, 1)
  D <- matrix(1, nrow = n_strain, ncol = 1)
  D[n_strain, 1] <- 0
  A <- matrix(1, nrow = 1, ncol = n_strain)

  model_obj <- list(
    y = y, N = 1L, P = 1L,
    llik_tab = snp.slicer:::build_categorical_llik_tab(0.05, 0.05)
  )
  state <- list(A = A, D = D, kmin = 1L, mixed = 1L)
  state$loglik <- snp.slicer:::categorical_loglikelihood_matrix(A, D, model_obj)

  prop <- (A %*% D) / rowSums(A)
  expect_gt(prop[1, 1], 0.99)      # in the offending window ...
  expect_lt(prop[1, 1], 1)         # ... but not equal to 1
  expect_true(is.infinite(state$loglik))

  fixed <- snp.slicer:::categorical_resolve_exceptions(state, model_obj)
  expect_false(is.infinite(fixed$loglik))
})

test_that("categorical_prop_bin_vec agrees with the scalar bin function", {
  probs <- c(-0.1, 0, 0.001, 0.5, 0.99, 0.991, 0.999, 1)
  expect_equal(
    snp.slicer:::categorical_prop_bin_vec(probs),
    vapply(probs, snp.slicer:::categorical_prop_bin, integer(1))
  )
})
