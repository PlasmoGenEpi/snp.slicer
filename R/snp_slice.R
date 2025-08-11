#' Bayesian Nonparametric Resolution of Multi-Strain Infections
#'
#' @description
#' SNP-Slice is a Bayesian nonparametric method for resolving multi-strain infections
#' using slice sampling with stick-breaking construction. The algorithm simultaneously
#' unveils strain haplotypes and links them to hosts from sequencing data.
#'
#' @param data Input data. Can be a matrix, data.frame, or file path. For read count data,
#'   should be a list with elements `read1` and `read0` (or `total`). For categorical data,
#'   should be a matrix with values 0, 0.5, or 1.
#' @param model Observation model to use. Options: "categorical", "poisson", "binomial",
#'   "negative_binomial" (default).
#' @param n_mcmc Number of MCMC iterations (default: 10000).
#' @param burnin Burn-in period. If NULL, defaults to n_mcmc/2.
#' @param alpha IBP concentration parameter (default: 2.6).
#' @param rho Dictionary sparsity parameter (default: 0.5).
#' @param threshold Threshold for identifying single infections (default: 0.001).
#' @param gap Early stopping threshold. If NULL, runs for full n_mcmc iterations.
#' @param seed Random seed for reproducibility.
#' @param verbose Whether to print progress information (default: TRUE).
#' @param store_mcmc Whether to store full MCMC samples (default: FALSE).
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
#' result <- snp_slice(data, model = "negative_binomial", n_mcmc = 1000)
#'
#' # Extract results
#' strains <- extract_strains(result)
#' allocations <- extract_allocations(result)
#' }
#'
#' @export
snp_slice <- function(data,
                      model = "negative_binomial",
                      n_mcmc = 10000,
                      burnin = NULL,
                      alpha = 2.6,
                      rho = 0.5,
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
  validate_input_data(data, model)
  validate_parameters(alpha, rho, threshold)
  validate_mcmc_settings(n_mcmc, burnin, gap)

  # Set default burnin if not provided
  if (is.null(burnin)) {
    burnin <- floor(n_mcmc / 2)
  }

  # Preprocess data
  processed_data <- preprocess_data(data, model)

  # Create model object
  model_obj <- create_model(model, processed_data, alpha = alpha, rho = rho, ...)

  # Run MCMC
  if (verbose) {
    cat("Running SNP-Slice with", model, "model\n")
    cat("N =", nrow(processed_data$y), "hosts, P =", ncol(processed_data$y), "SNPs\n")
    cat("MCMC iterations:", n_mcmc, "burnin:", burnin, "\n")
  }

  result <- run_mcmc(
    model_obj = model_obj,
    n_mcmc = n_mcmc,
    burnin = burnin,
    gap = gap,
    verbose = verbose,
    store_mcmc = store_mcmc
  )

  # Create results object
  results <- create_results_object(result, model_obj, processed_data)

  if (verbose) {
    cat("Analysis complete\n")
    if (log_performance) {
      cat("Performance Summary:\n")
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
