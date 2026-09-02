# Tests for calculate_allele_frequencies and calculate_allele_frequencies_by_sets

# Minimal snp_slice_results with processed_data required for allele frequency functions
make_allele_freq_result <- function() {
  A <- matrix(c(2, 0, 1, 1, 1, 0), nrow = 2, ncol = 3)  # 2 hosts, 3 strains
  D <- matrix(c(1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L), nrow = 3, ncol = 3)  # 3 strains, 3 SNPs
  P <- ncol(D)
  list(
    map_allocation_matrix = A,
    map_dictionary_matrix = D,
    # Default to a final sample identical to the MAP so that value-based tests
    # using the default estimate are unaffected; tests that need them to differ
    # override these fields.
    final_allocation_matrix = A,
    final_dictionary_matrix = D,
    model_info = list(
      model = "negative_binomial",
      N = 2,
      P = P,
      data_type = "read_counts",
      processed_data = list(
        target_ids = paste0("target_", seq_len(P)),
        r0_values = rep("ref", P),
        r1_values = rep("alt", P)
      )
    )
  ) |> structure(class = "snp_slice_results")
}

test_that("calculate_allele_frequencies returns correct structure", {
  result <- make_allele_freq_result()
  snp_indices <- c(1L, 2L, 3L)
  freq_df <- calculate_allele_frequencies(result, snp_indices)

  expect_s3_class(freq_df, "data.frame")
  expect_equal(
    sort(colnames(freq_df)),
    sort(c("allele", "frequency", "count", "total_parasites"))
  )
  expect_true(all(freq_df$frequency >= 0 & freq_df$frequency <= 1))
  expect_equal(sum(freq_df$frequency), 1)
  expect_true(all(freq_df$count >= 0))
  expect_equal(sum(freq_df$count), unique(freq_df$total_parasites))
  expect_true(nrow(freq_df) >= 1)
})

test_that("calculate_allele_frequencies with character indices matches integer indices when target_ids exist", {
  result <- make_allele_freq_result()
  tid <- result$model_info$processed_data$target_ids
  int_freq <- calculate_allele_frequencies(result, c(1L, 2L, 3L))
  char_freq <- calculate_allele_frequencies(result, c(tid[1], tid[2], tid[3]))

  expect_equal(int_freq$allele, char_freq$allele)
  expect_equal(int_freq$frequency, char_freq$frequency)
  expect_equal(int_freq$count, char_freq$count)
  expect_equal(int_freq$total_parasites, char_freq$total_parasites)
})

test_that("calculate_allele_frequencies errors on invalid inputs", {
  result <- make_allele_freq_result()
  expect_error(calculate_allele_frequencies(result, integer(0)), "non-empty")
  expect_error(calculate_allele_frequencies(result, c(1L, 1L)), "unique")
  ncol_d <- ncol(result$map_dictionary_matrix)
  expect_error(calculate_allele_frequencies(result, c(1L, ncol_d + 1L)), "valid SNP positions")
})

test_that("calculate_allele_frequencies with estimate = 'final_sample' uses the final matrices", {
  result <- make_allele_freq_result()
  # A final sample that differs from the MAP: shift a host onto another strain
  result$final_allocation_matrix <- matrix(c(2, 0, 0, 1, 1, 1), nrow = 2, ncol = 3)

  final_df <- calculate_allele_frequencies(result, c(1L, 2L, 3L), estimate = "final_sample")
  map_df <- calculate_allele_frequencies(result, c(1L, 2L, 3L), estimate = "map")

  # final sample is the default
  expect_equal(final_df, calculate_allele_frequencies(result, c(1L, 2L, 3L)))
  # same column shape as MAP
  expect_equal(sort(colnames(final_df)), sort(colnames(map_df)))
  expect_equal(sort(colnames(final_df)), sort(c("allele", "frequency", "count", "total_parasites")))
  # the two point estimates disagree because the allocations differ
  expect_false(isTRUE(all.equal(final_df, map_df)))
})

test_that("calculate_allele_frequencies with estimate = 'posterior' when MCMC samples exist", {
  result <- make_allele_freq_result()
  result$mcmc_samples <- list(
    list(A = result$map_allocation_matrix, D = result$map_dictionary_matrix),
    list(A = result$map_allocation_matrix, D = result$map_dictionary_matrix),
    list(A = result$map_allocation_matrix, D = result$map_dictionary_matrix),
    list(A = result$map_allocation_matrix, D = result$map_dictionary_matrix),
    list(A = result$map_allocation_matrix, D = result$map_dictionary_matrix)
  )
  freq_df <- calculate_allele_frequencies(result, c(1L, 2L), estimate = "posterior", n_samples = 5)
  expect_s3_class(freq_df, "data.frame")
  mcmc_cols <- c("allele", "frequency", "frequency_sd", "frequency_lower", "frequency_upper", "mean_count", "n_samples")
  expect_equal(sort(colnames(freq_df)), sort(mcmc_cols))
  expect_true(all(freq_df$frequency >= 0 & freq_df$frequency <= 1))
  expect_equal(sum(freq_df$frequency), 1)
  expect_equal(unique(freq_df$n_samples), 5)
  expect_true(all(freq_df$mean_count >= 0))
  expect_false(is.null(attr(freq_df, "mean_total_parasites")))
})

test_that("calculate_allele_frequencies errors on invalid interval", {
  result <- make_allele_freq_result()
  expect_error(calculate_allele_frequencies(result, c(1L, 2L), interval = 0), "interval must be")
  expect_error(calculate_allele_frequencies(result, c(1L, 2L), interval = 1), "interval must be")
})

test_that("calculate_allele_frequencies_by_sets returns list of frequency tables", {
  result <- make_allele_freq_result()
  target_sets <- list(single = c(1L), multi = c(1L, 2L, 3L))
  freqs <- calculate_allele_frequencies_by_sets(result, target_sets)

  expect_type(freqs, "list")
  expect_equal(length(freqs), 2)
  expect_equal(names(freqs), c("single", "multi"))

  for (i in seq_along(freqs)) {
    df <- freqs[[i]]
    expect_s3_class(df, "data.frame")
    expect_equal(sort(colnames(df)), sort(c("allele", "frequency", "count", "total_parasites")))
    expect_equal(sum(df$frequency), 1)
  }
})

test_that("calculate_allele_frequencies_by_sets uses set_1, set_2 when target_sets unnamed", {
  result <- make_allele_freq_result()
  target_sets <- list(c(1L), c(2L, 3L))
  freqs <- calculate_allele_frequencies_by_sets(result, target_sets)

  expect_equal(names(freqs), c("set_1", "set_2"))
  expect_equal(sum(freqs$set_1$frequency), 1)
  expect_equal(sum(freqs$set_2$frequency), 1)
})

test_that("calculate_allele_frequencies_by_sets errors on invalid target_sets", {
  result <- make_allele_freq_result()
  expect_error(calculate_allele_frequencies_by_sets(result, "not a list"), "must be a list")
  expect_error(calculate_allele_frequencies_by_sets(result, list()), "at least one set")
  expect_error(
    calculate_allele_frequencies_by_sets(result, list(integer(0))),
    "non-empty vector"
  )
})

# Regression: the dictionary-to-allele-label lookup is indexed by D + 1, and the
# two model families assign the slots in opposite directions. Count models put
# the minor allele in slot 1 and set D to 1 for it, so D == 0 denotes r1. The
# categorical model encodes slot 1 as observation 0, so D == 0 denotes r0.
make_slot_result <- function(model) {
  # One host carrying one strain; that strain is 0 at the single SNP.
  A <- matrix(1L, nrow = 1, ncol = 1)
  D <- matrix(0L, nrow = 1, ncol = 1)
  list(
    map_allocation_matrix = A,
    map_dictionary_matrix = D,
    final_allocation_matrix = A,
    final_dictionary_matrix = D,
    model_info = list(
      model = model,
      N = 1,
      P = 1,
      data_type = if (model == "categorical") "categorical" else "read_counts",
      processed_data = list(
        target_ids = "target_1",
        r0_values = "slot1",
        r1_values = "slot2"
      )
    )
  ) |> structure(class = "snp_slice_results")
}

test_that("a count-model dictionary entry of 0 is labelled with the slot-2 allele", {
  af <- calculate_allele_frequencies(make_slot_result("negative_binomial"), 1)
  top <- af[af$frequency > 0, ]
  expect_equal(nrow(top), 1L)
  expect_equal(top$allele, "slot2")
})

test_that("a categorical dictionary entry of 0 is labelled with the slot-1 allele", {
  af <- calculate_allele_frequencies(make_slot_result("categorical"), 1)
  top <- af[af$frequency > 0, ]
  expect_equal(nrow(top), 1L)
  expect_equal(top$allele, "slot1")
})

test_that("the two model families label the same dictionary in opposite directions", {
  count_af <- calculate_allele_frequencies(make_slot_result("poisson"), 1)
  cat_af <- calculate_allele_frequencies(make_slot_result("categorical"), 1)
  expect_false(identical(
    count_af$allele[count_af$frequency > 0],
    cat_af$allele[cat_af$frequency > 0]
  ))
})
