#' Bayesian Nonparametric Resolution of Multi-Strain Infections
#'
#' @description
#' SNP-Slice is a Bayesian nonparametric method for resolving multi-strain infections
#' using slice sampling with stick-breaking construction. The algorithm simultaneously
#' unveils strain haplotypes and links them to hosts from sequencing data.
#'
#' @param data Input data. Can be a matrix, data.frame, or file path. For read count data,
#'   should be a list with elements `read1` and `read0` (or `total`). For categorical data,
#'   can be a matrix with values 0, 0.5, or 1; or a long-format data.frame with columns
#'   \code{specimen_id}, \code{target_id}, \code{target_value}, and \code{target_count}.
#'   For a categorical data.frame, counts are converted to categories: ref-only -> 0,
#'   alt-only -> 1, both present -> 0.5, zero total -> NA. Matrix and categorical file
#'   inputs (e.g. \code{*_cat.txt}) remain supported.
#' @param model Observation model to use. Options: "categorical", "poisson", "binomial",
#'   "negative_binomial" (default).
#' @param n_sample Number of post-burn-in iterations to retain (default: 10000).
#'   Burn-in iterations are additional: the chain runs \code{n_burnin + n_sample}
#'   iterations in total and only the last \code{n_sample} are retained.
#' @param n_burnin Number of iterations discarded before sampling begins. If
#'   NULL, defaults to \code{floor(n_sample / 2)}.
#' @param alpha IBP concentration parameter (default: 2.6).
#' @param rho Dictionary sparsity parameter (default: 0.5 for categorical model and NULL otherwise, which means use the global minor allele frequency)
#' @param threshold Threshold for identifying single infections (default: 0.001).
#' @param gap Early stopping threshold. If NULL, runs for all
#'   \code{n_burnin + n_sample} iterations.
#' @param seed Random seed for reproducibility.
#' @param verbose Whether to print progress information (default: TRUE).
#' @param log_performance Whether to log performance metrics (default: FALSE).
#' @param store_mcmc Whether to store full MCMC samples (default: FALSE).
#'   Only post-burn-in iterations are stored.
#' @param ... Additional model-specific parameters.
#'
#' @return An object of class `snp_slice_results` containing:
#'   - `allocation_matrix`: Binary allocation matrix (A)
#'   - `dictionary_matrix`: Binary dictionary matrix (D)
#'   - `mcmc_samples`: MCMC samples (if store_mcmc = TRUE)
#'   - `diagnostics`: Convergence diagnostics
#'   - `parameters`: Model parameters used
#'   - `model_info`: Model specification
#'
#' @importFrom stats runif dpois dbinom dnbinom rbeta acf median
#' @importFrom utils read.delim tail
#'
#' @examples
#' \dontrun{
#' # Example with read count data
#' data <- list(
#'   read1 = matrix(c(10, 5, 15, 8), nrow = 2),
#'   read0 = matrix(c(90, 95, 85, 92), nrow = 2)
#' )
#'
#' result <- snp_slice(data, model = "negative_binomial", n_sample = 1000)
#'
#' # Extract results
#' strains <- extract_strains(result)
#' allocations <- extract_allocations(result)
#' }
#'
#' @export
snp_slice <- function(data,
                      model = "negative_binomial",
                      n_sample = 10000,
                      n_burnin = NULL,
                      alpha = 2.6,
                      rho = if (model == "categorical") 0.5 else NULL,
                      threshold = 0.001,
                      gap = NULL,
                      seed = NULL,
                      verbose = TRUE,
                      log_performance = FALSE,
                      store_mcmc = FALSE,
                      ...) {

  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  validate_parameters(alpha, rho, threshold)
  validate_mcmc_settings(n_sample, n_burnin, gap)

  # Set default burn-in if not provided
  if (is.null(n_burnin)) {
    n_burnin <- floor(n_sample / 2)
  }

  # Validate input data before preprocessing
  validate_input_data(data, model, ...)

  # Preprocess data
  processed_data <- preprocess_data(data, model, ...)

  # Create model object
  model_obj <- create_model(model, processed_data, alpha = alpha, rho = rho, ...)

  # Run MCMC
  if (verbose) {
    cat("Running SNP-Slice with", model, "model\n")
    cat("N =", nrow(processed_data$y), "hosts, P =", ncol(processed_data$y), "SNPs\n")
    cat("Retained samples:", n_sample, "burn-in:", n_burnin, "\n")
  }

  result <- run_mcmc(
    model_obj = model_obj,
    n_sample = n_sample,
    n_burnin = n_burnin,
    gap = gap,
    verbose = verbose,
    store_mcmc = store_mcmc
  )

  # Create results object
  results <- create_results_object(result, model_obj, processed_data)

  if (verbose) {
    cat("Analysis complete\n")
    if (log_performance) {
      print_performance_summary()
    }
  }

  return(results)
}

#' Model-specific constructors
#'
#' @rdname snp_slice
#' @param e1 Error parameter for categorical model (default: 0.05)
#' @param e2 Error parameter for categorical model (default: 0.05)
#' @export
snp_slice_categorical <- function(data, e1 = 0.05, e2 = 0.05, ...) {
  snp_slice(data, model = "categorical", e1 = e1, e2 = e2, ...)
}

#' @rdname snp_slice
#' @export
snp_slice_poisson <- function(data, ...) {
  snp_slice(data, model = "poisson", ...)
}

#' @rdname snp_slice
#' @export
snp_slice_binomial <- function(data, ...) {
  snp_slice(data, model = "binomial", ...)
}

#' @rdname snp_slice
#' @export
snp_slice_negative_binomial <- function(data, ...) {
  snp_slice(data, model = "negative_binomial", ...)
}
