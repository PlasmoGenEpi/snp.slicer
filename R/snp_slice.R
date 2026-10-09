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
#' @param model Observation model to use. One of \code{"categorical"},
#'   \code{"poisson"}, \code{"binomial"}, \code{"negative_binomial"}
#'   (default), or \code{"multinomial"}. Only the multinomial model keeps
#'   targets with more than two alleles; the others drop them. See
#'   \sQuote{The multinomial model} below.
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
#' @param n_chains Number of independent MCMC chains to run (default: 1). Each
#'   chain is seeded separately; the chain reaching the highest MAP log
#'   posterior supplies the top-level estimates and all chains are kept in
#'   \code{result$chains}.
#' @param n_cores Number of cores used to run chains simultaneously
#'   (default: 1, i.e. chains run sequentially). Capped at \code{n_chains}.
#' @param seed Random seed for reproducibility. Per-chain seeds are based on 
#'   this seed, so a full multi-chain run is reproducible.
#' @param verbose Whether to print progress information (default: TRUE).
#' @param log_performance Whether to log performance metrics (default: FALSE).
#' @param store_mcmc Whether to store full MCMC samples (default: FALSE).
#'   Only post-burn-in iterations are stored.
#' @param ... Additional model-specific parameters.
#'
#' @return An object of class `snp_slice_results` containing:
#'   - `chains`: Per-chain results, all stored the same way. Each holds that
#'     chain's MAP estimate (`map_allocation_matrix` (A), `map_dictionary_matrix`
#'     (D)) and final-sample estimate (`final_allocation_matrix`,
#'     `final_dictionary_matrix`), plus `mcmc_samples` (if store_mcmc = TRUE),
#'     `diagnostics`, `convergence`, and `seed`
#'   - `best_chain`: Index of the chain with the highest MAP log posterior
#'   - `parameters`: MCMC settings used
#'   - `model_info`: Model specification
#'
#'   The object holds no estimates of its own. Reach a chain's estimates with
#'   [get_chain()], which defaults to the best chain, or with
#'   [extract_allocations()] / [extract_strains()]; every diagnostic function
#'   also takes a `chain` argument. [compare_chains()] summarises all chains.
#'
#' @section The multinomial model:
#' The other models describe a target with a single number per specimen: the
#' dictionary is binary, so \code{(A \%*\% D) / rowSums(A)} gives the fraction of a
#' specimen's strains carrying the alternate allele, and the reference allele is
#' whatever is left over. That only works for two alleles, so targets with more
#' are dropped.
#'
#' The multinomial model instead describes a target with a probability
#' \emph{vector}, one entry per allele. The dictionary holds an allele code per
#' strain per target rather than a bit; expanding it to one indicator column per
#' (target, allele) pair makes each entry of
#' \code{(A \%*\% Dx) / rowSums(A)} the fraction of a specimen's strains carrying
#' that allele, and the entries for one target sum to one. The observed read
#' counts across a target's alleles are then multinomial in that vector. Nothing
#' limits the vector's length, so every allele at every target is kept.
#'
#' With exactly two alleles the vector is \code{(1 - p, p)} and the kernel reduces
#' to the binomial one, so the two models agree on biallelic data.
#'
#' Parameters read from \code{...}: \code{dict_prior} sets the per-target
#' categorical prior on dictionary entries, either \code{"empirical"} (default;
#' pooled allele read fractions with a pseudocount of one), \code{"uniform"}, or
#' a list of per-target probability vectors. \code{rho}, used by the count
#' models, does not apply.
#'
#' Implementation: the allocation update runs in compiled code as for the other
#' models, while the dictionary update is an R Gibbs step over alleles.
#'
#' @importFrom stats runif dpois dbinom dnbinom rbeta median
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
                      n_chains = 1,
                      n_cores = 1,
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
  validate_mcmc_settings(n_sample, n_burnin, gap, n_chains, n_cores)

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
    cat("Retained samples:", n_sample, "burn-in:", n_burnin,
        "chains:", n_chains, "\n")
  }

  result <- run_chains(
    model_obj = model_obj,
    n_sample = n_sample,
    n_burnin = n_burnin,
    gap = gap,
    verbose = verbose,
    store_mcmc = store_mcmc,
    n_chains = as.integer(n_chains),
    n_cores = as.integer(n_cores)
  )

  # Create results object
  results <- create_results_object(result, model_obj, processed_data)

  if (verbose) {
    cat("Analysis complete\n")
    if (log_performance) {
      # Chains may have run in separate processes, so use the timings the best
      # chain returned rather than this process's global log
      print_performance_summary(result$performance)
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

#' @rdname snp_slice
#' @param dict_prior Dictionary prior for the multinomial model; see
#'   \code{model}.
#' @export
snp_slice_multinomial <- function(data, dict_prior = "empirical", ...) {
  snp_slice(data, model = "multinomial", dict_prior = dict_prior, ...)
}
