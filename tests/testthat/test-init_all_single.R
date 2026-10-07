# Regression: the count-model initializers reorder strains so those belonging
# to single infections come first. When every strain already belongs to a
# single infection there is nothing to append, and the unguarded
# `ord[(n + 1):nstrain] <- setdiff(...)` failed with "replacement has length
# zero". That happens when all specimens are single infections, and also when
# each mixed specimen's rounded genotype coincides with some single specimen's.

all_single <- list(
  read0 = matrix(c(100, 0, 100, 0), 2, 2, dimnames = list(c("s1", "s2"), c("t1", "t2"))),
  read1 = matrix(c(0, 100, 0, 100), 2, 2, dimnames = list(c("s1", "s2"), c("t1", "t2")))
)

# s3 is mixed (50/50) but rounds to s1's genotype, so no strain is "new".
mixed_rounds_to_single <- list(
  read0 = matrix(c(100, 0, 50, 100, 0, 50), 3, 2, dimnames = list(c("s1", "s2", "s3"), c("t1", "t2"))),
  read1 = matrix(c(0, 100, 50, 0, 100, 50), 3, 2, dimnames = list(c("s1", "s2", "s3"), c("t1", "t2")))
)

for (model in c("binomial", "negative_binomial", "poisson")) {
  test_that(paste(model, "initializes when every strain belongs to a single infection"), {
    for (dat in list(all_single, mixed_rounds_to_single)) {
      pd <- snp.slicer:::preprocess_data(dat, model)
      m <- snp.slicer:::create_model(model, pd, alpha = 2.6, rho = 0.5)
      st <- m$initialize_state(m, 0.001)
      expect_true(is.finite(st$loglik))
      expect_equal(ncol(st$A), nrow(st$D))
      expect_false(anyNA(st$D))
      # Every specimen carries at least one strain.
      expect_true(all(rowSums(st$A) >= 1))
    }
    res <- snp_slice(all_single, model = model, n_sample = 10, n_burnin = 5, rho = 0.5,
                     seed = 1, verbose = FALSE)
    expect_s3_class(res, "snp_slice_results")
  })
}

test_that("categorical initializes with only single infections", {
  y <- matrix(c(0, 1, 0, 1), 2, 2, dimnames = list(c("s1", "s2"), c("t1", "t2")))
  res <- snp_slice(y, model = "categorical", n_sample = 10, n_burnin = 5, seed = 1, verbose = FALSE)
  expect_s3_class(res, "snp_slice_results")
})
