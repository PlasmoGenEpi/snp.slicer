# Tests for monomorphic-locus handling, long-format validation, and the
# duplicated-allele-label path through calculate_allele_counts_single.

test_that("load_dataframe duplicates the allele label for monomorphic loci (no '?')", {
  # t1 biallelic (A/T), t2 monomorphic (G only).
  df <- data.frame(
    specimen_id  = c("s1", "s1", "s1", "s2", "s2"),
    target_id    = c("t1", "t1", "t2", "t1", "t2"),
    target_value = c("A",  "T",  "G",  "A",  "G"),
    target_count = c(5,    3,    9,    4,    7)
  )
  out <- snp.slicer:::load_dataframe(df, model = "negative_binomial")

  expect_setequal(out$target_ids, c("t1", "t2"))
  t2 <- match("t2", out$target_ids)
  t1 <- match("t1", out$target_ids)

  # Monomorphic locus: both alleles labelled with the single observed value.
  expect_equal(out$r0_values[t2], "G")
  expect_equal(out$r1_values[t2], "G")

  # Biallelic locus is oriented by frequency: slot 1 (read0 / y) is the MINOR
  # allele. Here T (total 3) is the minor allele and A (total 9) the major.
  expect_equal(out$r0_values[t1], "T")
  expect_equal(out$r1_values[t1], "A")

  # No placeholder leaks into the labels.
  expect_false(any(out$r0_values == "?"))
  expect_false(any(out$r1_values == "?"))

  # Monomorphic locus: all reads land in read0 (y); read1 is a structural zero.
  # r is the total coverage (read0 + read1), so r == y == 9 confirms read1 == 0.
  expect_equal(out$y["s1", "t2"], 9)
  expect_equal(out$r["s1", "t2"], 9)
  expect_equal(out$y["s2", "t2"], 7)
  expect_equal(out$r["s2", "t2"], 7)
})

test_that("load_dataframe distinguishes missing genotypes from homozygous calls", {
  # t1 biallelic (A/T), t2 monomorphic (G). s2 has no t1 reads; s3 has no t2 reads.
  df <- data.frame(
    specimen_id  = c("s1", "s1", "s1", "s2", "s3", "s3"),
    target_id    = c("t1", "t1", "t2", "t2", "t1", "t1"),
    target_value = c("A",  "T",  "G",  "G",  "A",  "T"),
    target_count = c(5,    3,    9,    7,    2,    6)
  )
  out <- snp.slicer:::load_dataframe(df, model = "negative_binomial")

  # Missing genotypes -> NA (no reads at that locus for that specimen).
  expect_true(is.na(out$y["s3", "t2"]))
  expect_true(is.na(out$r["s3", "t2"]))
  expect_true(is.na(out$y["s2", "t1"]))

  # Homozygous monomorphic call is a real observation (read0 = 9, read1 = 0), not NA.
  expect_equal(out$y["s1", "t2"], 9)
  expect_equal(out$r["s1", "t2"], 9)
  expect_false(is.na(out$y["s1", "t2"]))
})

test_that("validate_input_data rejects duplicate (specimen, target, value) rows", {
  df <- data.frame(
    specimen_id  = c("s1", "s1", "s1", "s2"),
    target_id    = c("t1", "t1", "t1", "t1"),
    target_value = c("A",  "A",  "T",  "T"),   # (s1, t1, A) duplicated
    target_count = c(5,    2,    3,    4)
  )
  expect_error(
    snp.slicer:::validate_input_data(df, "negative_binomial"),
    "duplicate"
  )
})

test_that("validate_input_data rejects inputs with no biallelic locus", {
  df <- data.frame(
    specimen_id  = c("s1", "s2", "s3"),
    target_id    = c("t1", "t1", "t2"),
    target_value = c("A",  "A",  "G"),   # t1 and t2 both monomorphic
    target_count = c(5,    4,    3)
  )
  expect_error(
    snp.slicer:::validate_input_data(df, "negative_binomial"),
    "biallelic"
  )
})

test_that("calculate_allele_counts_single does not double-count duplicated allele labels", {
  # Single monomorphic SNP: r0 == r1 == "G".
  # ind1 carries strain1; ind2 carries strain1 + strain2. All strains are "G".
  A <- matrix(c(1, 0,
                1, 1), nrow = 2, byrow = TRUE)
  D <- matrix(c(0, 1), nrow = 2)  # strain1 bit 0, strain2 bit 1, one SNP
  res <- snp.slicer:::calculate_allele_counts_single(
    A, D, snp_indices = 1, r0_values = "G", r1_values = "G"
  )

  expect_equal(nrow(res), 1)
  expect_equal(res$allele, "G")
  expect_equal(res$count[1], 3)  # 1 + 2 strain-copies, all "G"
})

test_that("snp_slice runs end-to-end with a monomorphic locus present", {
  df <- data.frame(
    specimen_id  = c("s1", "s1", "s1", "s2", "s2", "s2"),
    target_id    = c("t1", "t1", "t2", "t1", "t1", "t2"),
    target_value = c("A",  "T",  "G",  "A",  "T",  "G"),
    target_count = c(5,    3,    9,    4,    6,    7)
  )
  result <- snp_slice(df, model = "negative_binomial", n_mcmc = 25, verbose = FALSE)
  expect_s3_class(result, "snp_slice_results")
})
