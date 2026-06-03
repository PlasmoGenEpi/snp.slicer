make_kplus2_mu_state <- function(model_obj) {
  state <- model_obj$initialize_state(model_obj, threshold = 0.001)
  state <- slice_init(state, model_obj)
  state$mu <- c(0.85, 0.25, 0.12, 0.06)
  state$kplus <- 2L
  state$ktrunc <- length(state$mu)
  state$kstar <- 1L
  state
}

test_that("slice_update_mu preserves mu[1] when kplus is 2", {
  processed <- preprocess_data(
    list(
      read1 = matrix(c(10, 5, 15, 8), nrow = 2),
      read0 = matrix(c(90, 95, 85, 92), nrow = 2)
    ),
    "negative_binomial"
  )
  model_obj <- create_model("negative_binomial", processed)
  state <- make_kplus2_mu_state(model_obj)

  set.seed(123)
  updated <- slice_update_mu(state, model_obj)

  expect_gt(updated$mu[1], updated$mu[2])
  expect_lt(updated$mu[1], 1)
  expect_false(abs(updated$mu[1] - 0.5) < 1e-12)
})
