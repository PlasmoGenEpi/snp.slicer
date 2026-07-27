# Global variables for data.frame column names
utils::globalVariables(c("iteration", "logpost", "kstar", "n_strains"))

#' Remove additional burn-in samples from MCMC results
#'
#' Stored MCMC samples already exclude burn-in, so this returns the full stored
#' chain unless \code{additional_burnin} asks for further trimming.
#'
#' @param results SNP-Slice results object
#' @param additional_burnin Number of additional stored samples to discard from
#'   the start of the retained chain
#'
#' @return List of MCMC samples
#' @keywords internal
retained_samples <- function(results, additional_burnin = 0) {
  if (is.null(results$mcmc_samples)) {
    stop("MCMC samples not stored. Set store_mcmc = TRUE when running snp_slice()")
  }
  if (!is.numeric(additional_burnin) || length(additional_burnin) != 1 ||
        additional_burnin < 0 || additional_burnin != as.integer(additional_burnin)) {
    stop("additional_burnin must be a non-negative integer")
  }

  n <- length(results$mcmc_samples)
  if (additional_burnin >= n) {
    stop("additional_burnin (", additional_burnin, ") leaves no samples; only ",
         n, " retained samples are available")
  }

  results$mcmc_samples[seq(additional_burnin + 1, n)]
}

#' Resolve the chain a diagnostic function should use
#'
#' Diagnostics take a \code{chain} argument rather than requiring a results
#' object that holds a single chain. \code{NULL} (the default everywhere) means
#' the best chain.
#'
#' @param results SNP-Slice results object
#' @param chain Chain index, or \code{NULL} for the best chain
#'
#' @return A single-chain results object, as returned by \code{\link{get_chain}}
#' @keywords internal
resolve_chain <- function(results, chain = NULL) {
  if (is.null(chain)) {
    chain <- results$best_chain
  }
  get_chain(results, chain)
}

#' Extract a single chain as a results object
#'
#' @description
#' A results object stores every chain identically in \code{results$chains} and
#' keeps no estimates of its own. This function flattens one chain into a
#' standalone \code{snp_slice_results} object carrying that chain's
#' \code{map_allocation_matrix}, \code{map_dictionary_matrix},
#' \code{final_allocation_matrix}, \code{final_dictionary_matrix}, and
#' \code{mcmc_samples}, which is how the estimates of a run are reached.
#' Diagnostic functions call it for you via their \code{chain} argument.
#'
#' @param results SNP-Slice results object
#' @param chain Index of the chain to extract (defaults to the best chain)
#'
#' @return An \code{snp_slice_results} object for the requested chain
#' @export
#' @examples
#' result <- load_example_results()
#' chain1 <- get_chain(result, 1)
#' summary(chain1)
#' dim(chain1$map_allocation_matrix)
get_chain <- function(results, chain = results$best_chain) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }

  if (is.null(results$chains)) {
    if (!is.null(chain) && !identical(as.integer(chain), 1L)) {
      stop("chain must be 1: this object already holds a single chain")
    }
    return(results)
  }

  chains <- results$chains
  if (is.null(chain)) {
    chain <- results$best_chain
  }
  if (!is.numeric(chain) || length(chain) != 1 || chain < 1 ||
        chain > length(chains)) {
    stop("chain must be an integer between 1 and ", length(chains))
  }

  chain_result <- chains[[chain]]

  out <- list(
    chain_id = chain_result$chain_id,
    seed = chain_result$seed,
    map_allocation_matrix = chain_result$map_allocation_matrix,
    map_dictionary_matrix = chain_result$map_dictionary_matrix,
    final_allocation_matrix = chain_result$final_allocation_matrix,
    final_dictionary_matrix = chain_result$final_dictionary_matrix,
    mcmc_samples = chain_result$mcmc_samples,
    diagnostics = chain_result$diagnostics,
    convergence = chain_result$convergence,
    performance = chain_result$performance,
    parameters = results$parameters,
    model_info = results$model_info
  )

  class(out) <- "snp_slice_results"
  out
}

#' Compare chains from a multi-chain run
#'
#' @param results SNP-Slice results object
#'
#' @return A data frame with one row per chain: chain, seed, iterations_run,
#'   map_logpost, final_logpost, map_kstar, n_strains, gap_converged, and
#'   whether the chain was selected as best
#' @export
#' @examples
#' result <- load_example_results()
#' compare_chains(result)
compare_chains <- function(results) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  if (is.null(results$chains)) {
    stop("results has no chains; it was produced by a version of snp.slicer ",
         "that predates multi-chain support. Re-run snp_slice() to use this function")
  }
  chains <- results$chains

  data.frame(
    chain = vapply(chains, function(x) as.integer(x$chain_id), integer(1)),
    seed = vapply(chains, function(x) as.numeric(x$seed), numeric(1)),
    iterations_run = vapply(chains,
                            function(x) as.integer(x$convergence$iterations_run), integer(1)),
    map_logpost = vapply(chains, function(x) x$diagnostics$map_logpost, numeric(1)),
    final_logpost = vapply(chains, function(x) x$diagnostics$final_logpost, numeric(1)),
    map_kstar = vapply(chains, function(x) as.integer(x$diagnostics$map_kstar), integer(1)),
    n_strains = vapply(chains, function(x) nrow(x$map_dictionary_matrix), integer(1)),
    gap_converged = vapply(chains,
                           function(x) isTRUE(x$convergence$gap_converged), logical(1)),
    best = seq_along(chains) == results$best_chain,
    stringsAsFactors = FALSE
  )
}

#' Extract strain information from SNP-Slice results
#'
#' @param results SNP-Slice results object
#' @inheritParams calculate_allele_frequencies
#'
#' @return List containing strain information
#' @export
extract_strains <- function(results, chain = NULL) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  results <- resolve_chain(results, chain)

  # Extract strain information
  strains <- list(
    dictionary = results$map_dictionary_matrix,
    n_strains = nrow(results$map_dictionary_matrix),
    n_snps = ncol(results$map_dictionary_matrix),
    strain_names = paste0("Strain_", 1:nrow(results$map_dictionary_matrix))
  )
  
  return(strains)
}

#' Extract allocation information from SNP-Slice results
#'
#' @param results SNP-Slice results object
#' @inheritParams calculate_allele_frequencies
#'
#' @return List containing allocation information
#' @export
extract_allocations <- function(results, chain = NULL) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  results <- resolve_chain(results, chain)

  # Extract allocation information (MAP estimate)
  n_hosts <- nrow(results$map_allocation_matrix)
  n_strains <- ncol(results$map_allocation_matrix)
  allocations <- list(
    allocation_matrix = results$map_allocation_matrix,
    n_hosts = n_hosts,
    n_strains = n_strains,
    host_names = paste0("Host_", seq_len(n_hosts)),
    strain_names = paste0("Strain_", seq_len(n_strains)),
    multiplicity_of_infection = rowSums(results$map_allocation_matrix)
  )
  
  return(allocations)
}

#' Calculate estimated individual COI with uncertainty
#'
#' @description
#' Returns per-host complexity of infection (COI). A point estimate
#' (\code{estimate = "final_sample"} or \code{"map"}) gives one COI per host with
#' \code{NA} uncertainty columns; \code{estimate = "posterior"} computes
#' posterior mean, SD, and credible interval from the MCMC samples.
#'
#' @param results A \code{snp_slice_results} object.
#' @inheritParams calculate_allele_frequencies
#'
#' @return A data frame with one row per host: host_index, host_id, coi_estimate,
#'   coi_sd, coi_lower, coi_upper. Uncertainty columns are NA for a point
#'   estimate or when no MCMC samples are available.
#'
#' @export
#' @examples
#' result <- load_example_results()
#' coi_final <- calculate_individual_coi(result, estimate = "final_sample")
#' head(coi_final)
#' if (!is.null(get_chain(result)$mcmc_samples)) {
#'   coi_post <- calculate_individual_coi(result, estimate = "posterior", n_samples = 50)
#'   head(coi_post)
#' }
calculate_individual_coi <- function(results,
                                    estimate = c("final_sample", "map", "posterior"),
                                    n_samples = 100,
                                    interval = 0.95,
                                    additional_burnin = 0,
                                    chain = NULL) {
  if (!inherits(results, "snp_slice_results")) {
    stop("results must be an snp_slice_results object")
  }
  estimate <- match.arg(estimate)
  results <- resolve_chain(results, chain)
  if (!is.numeric(n_samples) || n_samples < 1) {
    stop("n_samples must be a positive integer")
  }
  if (!is.numeric(interval) || interval <= 0 || interval >= 1) {
    stop("interval must be a number between 0 and 1 (exclusive)")
  }

  # Host count is invariant across estimates; use the always-present MAP matrix.
  n_hosts <- nrow(results$map_allocation_matrix)

  host_id <- NA_character_
  if (!is.null(results$model_info$processed_data$specimen_ids)) {
    sid <- results$model_info$processed_data$specimen_ids
    if (length(sid) == n_hosts) {
      host_id <- sid
    }
  }
  if (length(host_id) == 1L && is.na(host_id)) {
    host_id <- rep(NA_character_, n_hosts)
  }

  if (estimate != "posterior") {
    A <- point_estimate_matrices(results, estimate)$A
    coi_estimate <- rowSums(A)
    out <- data.frame(
      host_index = seq_len(n_hosts),
      host_id = host_id,
      coi_estimate = as.integer(coi_estimate),
      coi_sd = NA_real_,
      coi_lower = NA_real_,
      coi_upper = NA_real_,
      stringsAsFactors = FALSE
    )
    return(out)
  }

  if (is.null(results$mcmc_samples)) {
    stop("MCMC samples not available. Use estimate = \"map\" or \"final_sample\", or run snp_slice with store_mcmc = TRUE")
  }

  samples <- retained_samples(results, additional_burnin)
  n_avail <- length(samples)

  if (n_samples > n_avail) {
    warning("Requested ", n_samples, " samples but only ", n_avail, " retained. Using all retained samples.")
    n_samples <- n_avail
  }

  sample_indices <- sample.int(n_avail, n_samples, replace = FALSE)

  coi_matrix <- vapply(sample_indices, function(i) {
    rowSums(samples[[i]]$A)
  }, FUN.VALUE = numeric(n_hosts))

  probs <- c((1 - interval) / 2, 1 - (1 - interval) / 2)
  coi_estimate <- rowMeans(coi_matrix)
  coi_sd <- apply(coi_matrix, 1L, stats::sd)
  coi_lower <- apply(coi_matrix, 1L, stats::quantile, probs = probs[1], names = FALSE)
  coi_upper <- apply(coi_matrix, 1L, stats::quantile, probs = probs[2], names = FALSE)

  data.frame(
    host_index = seq_len(n_hosts),
    host_id = host_id,
    coi_estimate = coi_estimate,
    coi_sd = coi_sd,
    coi_lower = coi_lower,
    coi_upper = coi_upper,
    stringsAsFactors = FALSE
  )
}

#' Plot convergence diagnostics
#'
#' @param results SNP-Slice results object
#' @param type Type of plot ("logpost", "kstar", "n_strains")
#' @inheritParams calculate_allele_frequencies
#'
#' @return Plot object (if ggplot2 is available)
#' @export
plot_convergence <- function(results, type = "logpost", additional_burnin = 0, chain = NULL) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  results <- resolve_chain(results, chain)

  if (!type %in% c("logpost", "kstar", "n_strains")) {
    stop("Invalid plot type. Choose from: 'logpost', 'kstar', 'n_strains'")
  }

  # Extract retained (post-burn-in) MCMC samples
  samples <- retained_samples(results, additional_burnin)

  # True chain positions, so the x-axis starts after burn-in
  iterations <- vapply(samples, function(s) s$iteration, numeric(1))

  if (type == "logpost") {
    logpost_values <- sapply(samples, function(s) s$logpost)
    plot_data <- data.frame(
      iteration = iterations,
      logpost = logpost_values
    )
    
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration, y = logpost)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "Log Posterior Convergence",
                     x = "MCMC Iteration",
                     y = "Log Posterior") +
        ggplot2::theme_minimal()
      return(p)
    } else {
      plot(plot_data$iteration, plot_data$logpost, type = "l",
           main = "Log Posterior Convergence",
           xlab = "MCMC Iteration", ylab = "Log Posterior")
    }
  } else if (type == "kstar") {
    kstar_values <- sapply(samples, function(s) s$kstar)
    plot_data <- data.frame(
      iteration = iterations,
      kstar = kstar_values
    )
    
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration, y = kstar)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "Number of Active Strains",
                     x = "MCMC Iteration",
                     y = "k* (Active Strains)") +
        ggplot2::theme_minimal()
      return(p)
    } else {
      plot(plot_data$iteration, plot_data$kstar, type = "l",
           main = "Number of Active Strains",
           xlab = "MCMC Iteration", ylab = "k* (Active Strains)")
    }
  } else if (type == "n_strains") {
    n_strains_values <- sapply(samples, function(s) sum(colSums(s$A) > 0))
    plot_data <- data.frame(
      iteration = iterations,
      n_strains = n_strains_values
    )
    
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration, y = n_strains)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "Number of Strains",
                     x = "MCMC Iteration",
                     y = "Number of Strains") +
        ggplot2::theme_minimal()
      return(p)
    } else {
      plot(plot_data$iteration, plot_data$n_strains, type = "l",
           main = "Number of Strains",
           xlab = "MCMC Iteration", ylab = "Number of Strains")
    }
  }
}

#' Print summary of SNP-Slice results
#'
#' @param object SNP-Slice results object
#' @inheritParams calculate_allele_frequencies
#' @param ... Additional arguments
#'
#' @return Summary information
#' @export
#' @importFrom stats median
summary.snp_slice_results <- function(object, chain = NULL, ...) {
  if (!inherits(object, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  # Estimates come from one chain; the chain overview needs the whole object
  results <- object
  object <- resolve_chain(results, chain)

  cat("SNP-Slice Results Summary\n")
  cat("========================\n\n")
  
  # Model information
  pd <- object$model_info$processed_data
  N <- if (!is.null(pd$N)) pd$N else object$model_info$N
  P <- if (!is.null(pd$P)) pd$P else object$model_info$P
  cat("Model:", object$model_info$model, "\n")
  cat("Data dimensions:", N, "hosts x", P, "SNPs\n")
  data_type <- object$model_info$data_type
  if (is.null(data_type)) data_type <- object$model_info$model
  cat("Data type:", data_type, "\n\n")
  
  # Results summary
  cat("Results:\n")
  cat("- Number of strains identified:", nrow(object$map_dictionary_matrix), "\n")
  cat("- Number of hosts:", nrow(object$map_allocation_matrix), "\n")

  # Multiplicity of infection
  moi <- rowSums(object$map_allocation_matrix)
  cat("- Multiplicity of infection (MOI):\n")
  cat("  - Mean MOI:", round(mean(moi), 2), "\n")
  cat("  - Median MOI:", round(stats::median(moi), 2), "\n")
  cat("  - Range:", min(moi), "-", max(moi), "\n")
  
  # Single vs mixed infections
  single_infections <- sum(moi == 1)
  mixed_infections <- sum(moi > 1)
  cat("  - Single infections:", single_infections, "(", round(100 * single_infections / length(moi), 1), "%)\n")
  cat("  - Mixed infections:", mixed_infections, "(", round(100 * mixed_infections / length(moi), 1), "%)\n\n")
  
  # Chain information
  n_chains <- length(results$chains)
  if (n_chains > 1) {
    cat("Chains:", n_chains, "(best chain:", results$best_chain, ")\n")
    print(compare_chains(results), row.names = FALSE)
    cat("\n")
  }

  # Convergence information
  if (!is.null(object$convergence)) {
    cat("Convergence:\n")
    cat("- Iterations run:", object$convergence$iterations_run, "\n")
    if (!is.null(object$convergence$samples_retained)) {
      cat("- Samples retained (post-burn-in):", object$convergence$samples_retained, "\n")
    }
    cat("- Gap Converged:", ifelse(object$convergence$gap_converged, "Yes", "No"), "\n")
  }
  
  # Diagnostics
  if (!is.null(object$diagnostics)) {
    cat("- Final log posterior:", round(object$diagnostics$final_logpost, 2), "\n")
    cat("- MAP log posterior:", round(object$diagnostics$map_logpost, 2), "\n")
    cat("- Final k*:", object$diagnostics$final_kstar, "\n")
    cat("- MAP k*:", object$diagnostics$map_kstar, "\n")
  }
  
  cat("\n")
}

#' Print SNP-Slice results
#'
#' @param x SNP-Slice results object
#' @inheritParams calculate_allele_frequencies
#' @param ... Additional arguments
#'
#' @return Print information
#' @export
print.snp_slice_results <- function(x, chain = NULL, ...) {
  if (!inherits(x, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  n_chains <- length(x$chains)
  best_chain <- x$best_chain
  x <- resolve_chain(x, chain)

  cat("SNP-Slice Results\n")
  cat("================\n")
  cat("Model:", x$model_info$model, "\n")
  cat("Dimensions:", nrow(x$map_allocation_matrix), "hosts x", ncol(x$map_dictionary_matrix), "strains x", ncol(x$map_dictionary_matrix), "SNPs\n")
  cat("Strains identified:", nrow(x$map_dictionary_matrix), "\n")

  if (n_chains > 1) {
    cat("Chains:", n_chains, "(best chain:", best_chain, ")\n")
  }

  if (!is.null(x$convergence)) {
    cat("Gap Converged:", ifelse(x$convergence$gap_converged, "Yes", "No"), "\n")
  }
  
  cat("\n")
}

