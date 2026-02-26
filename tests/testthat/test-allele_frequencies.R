# Tests for calculate_allele_frequencies and calculate_allele_frequencies_by_sets

# Minimal snp_slice_results with processed_data required for allele frequency functions
make_allele_freq_result <- function() {
  A <- matrix(c(2, 0, 1, 1, 1, 0), nrow = 2, ncol = 3)  # 2 hosts, 3 strains
  D <- matrix(c(1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L), nrow = 3, ncol = 3)  # 3 strains, 3 SNPs
  P <- ncol(D)
  list(
    allocation_matrix = A,
    dictionary_matrix = D,
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
  ncol_d <- ncol(result$dictionary_matrix)
  expect_error(calculate_allele_frequencies(result, c(1L, ncol_d + 1L)), "valid SNP positions")
})

test_that("calculate_allele_frequencies with use_map = FALSE when MCMC samples exist", {
  result <- make_allele_freq_result()
  result$mcmc_samples <- list(
    list(A = result$allocation_matrix, D = result$dictionary_matrix),
    list(A = result$allocation_matrix, D = result$dictionary_matrix),
    list(A = result$allocation_matrix, D = result$dictionary_matrix),
    list(A = result$allocation_matrix, D = result$dictionary_matrix),
    list(A = result$allocation_matrix, D = result$dictionary_matrix)
  )
  freq_df <- calculate_allele_frequencies(result, c(1L, 2L), use_map = FALSE, n_samples = 5)
  expect_s3_class(freq_df, "data.frame")
  expect_equal(sort(colnames(freq_df)), sort(c("allele", "frequency", "count", "total_parasites")))
  expect_true(all(freq_df$frequency >= 0 & freq_df$frequency <= 1))
  expect_equal(sum(freq_df$frequency), 1)
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
