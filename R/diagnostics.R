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

#' Extract strain information from SNP-Slice results
#'
#' @param results SNP-Slice results object
#'
#' @return List containing strain information
#' @export
extract_strains <- function(results) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  
  # Extract strain information
  strains <- list(
    dictionary = results$dictionary_matrix,
    n_strains = nrow(results$dictionary_matrix),
    n_snps = ncol(results$dictionary_matrix),
    strain_names = paste0("Strain_", 1:nrow(results$dictionary_matrix))
  )
  
  return(strains)
}

#' Extract allocation information from SNP-Slice results
#'
#' @param results SNP-Slice results object
#'
#' @return List containing allocation information
#' @export
extract_allocations <- function(results) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  
  # Extract allocation information
  n_hosts <- nrow(results$allocation_matrix)
  n_strains <- ncol(results$allocation_matrix)
  allocations <- list(
    allocation_matrix = results$allocation_matrix,
    n_hosts = n_hosts,
    n_strains = n_strains,
    host_names = paste0("Host_", seq_len(n_hosts)),
    strain_names = paste0("Strain_", seq_len(n_strains)),
    multiplicity_of_infection = rowSums(results$allocation_matrix)
  )
  
  return(allocations)
}

#' Calculate estimated individual COI with uncertainty
#'
#' @description
#' Returns per-host complexity of infection (COI). With \code{use_map = TRUE}
#' you get a point estimate (MAP); with \code{use_map = FALSE} and MCMC samples,
#' posterior mean, SD, and credible interval are computed.
#'
#' @param results A \code{snp_slice_results} object.
#' @param use_map If \code{TRUE} (default), use MAP only; uncertainty columns
#'   are \code{NA}. If \code{FALSE}, use MCMC samples for mean, SD, and interval.
#' @param n_samples When \code{use_map = FALSE}, number of MCMC samples to use
#'   (capped at the number of retained samples).
#' @param interval Numeric in (0, 1). Credible interval width when using MCMC
#'   (e.g. 0.95 for 2.5 and 97.5 percent quantiles).
#' @inheritParams calculate_allele_frequencies
#'
#' @return A data frame with one row per host: host_index, host_id, coi_estimate,
#'   coi_sd, coi_lower, coi_upper. Uncertainty columns are NA when using MAP or
#'   when no MCMC samples are available.
#'
#' @export
#' @examples
#' result <- load_example_results()
#' coi_map <- calculate_individual_coi(result, use_map = TRUE)
#' head(coi_map)
#' if (!is.null(result$mcmc_samples)) {
#'   coi_post <- calculate_individual_coi(result, use_map = FALSE, n_samples = 50)
#'   head(coi_post)
#' }
calculate_individual_coi <- function(results,
                                    use_map = TRUE,
                                    n_samples = 100,
                                    interval = 0.95,
                                    additional_burnin = 0) {
  if (!inherits(results, "snp_slice_results")) {
    stop("results must be an snp_slice_results object")
  }
  if (!is.logical(use_map) || length(use_map) != 1) {
    stop("use_map must be a single logical value")
  }
  if (!is.numeric(n_samples) || n_samples < 1) {
    stop("n_samples must be a positive integer")
  }
  if (!is.numeric(interval) || interval <= 0 || interval >= 1) {
    stop("interval must be a number between 0 and 1 (exclusive)")
  }

  A <- results$allocation_matrix
  n_hosts <- nrow(A)

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

  if (use_map) {
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
    stop("MCMC samples not available. Set use_map = TRUE or run snp_slice with store_mcmc = TRUE")
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
plot_convergence <- function(results, type = "logpost", additional_burnin = 0) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }

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

#' Calculate effective sample size for MCMC diagnostics
#'
#' @param results SNP-Slice results object
#' @param parameter Parameter to analyze. Options include:
#'   - "logpost": Log posterior probability
#'   - "kstar": Number of active strains
#'   - "n_strains": Number of strains with non-zero allocations
#'   - "ktrunc": Truncation level
#'   - "mu": Stick-breaking weights (returns ESS for each weight)
#'   - "A": Allocation matrix (returns overall ESS)
#'   - "D": Dictionary matrix (returns overall ESS)
#'   - "all": Calculate ESS for all available parameters
#' @param method Method for ESS calculation. Options:
#'   - "autocorrelation": Standard autocorrelation-based ESS (default)
#'   - "batch_means": Batch means method
#'   - "spectral": Spectral density method
#' @inheritParams calculate_allele_frequencies
#'
#' @return Effective sample size(s) and diagnostic information
#' @export
#' @importFrom stats acf var
effective_sample_size <- function(results, parameter = "logpost", method = "autocorrelation",
                                  additional_burnin = 0) {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }

  samples <- retained_samples(results, additional_burnin)
  n_samples <- length(samples)

  if (n_samples < 2) {
    warning("Insufficient samples for ESS calculation")
    return(list(ess = n_samples, n_samples = n_samples, method = method))
  }
  
  # Helper function to calculate ESS using different methods
  calculate_ess <- function(values, method) {
    n <- length(values)
    if (n < 2) return(list(ess = n, acf_sum = 0, var_ratio = 1))

    if (method == "autocorrelation") {
      # Calculate autocorrelation
      max_lag <- min(n - 1, floor(n/4))  # Use at most 25% of sample size
      acf_result <- stats::acf(values, lag.max = max_lag, plot = FALSE)
      acf_values <- acf_result$acf[, 1, 1]
      
      # Sum of autocorrelations (excluding lag 0)
      acf_sum <- sum(acf_values[-1])
      
      # Calculate ESS
      ess <- n / (1 + 2 * acf_sum)
      
      return(list(ess = ess, acf_sum = acf_sum, var_ratio = 1 + 2 * acf_sum))
      
    } else if (method == "batch_means") {
      # Batch means method
      batch_size <- max(1, floor(sqrt(n)))
      n_batches <- floor(n / batch_size)
      
      if (n_batches < 2) {
        return(list(ess = n, acf_sum = 0, var_ratio = 1))
      }
      
      # Calculate batch means
      batch_means <- numeric(n_batches)
      for (i in 1:n_batches) {
        start_idx <- (i - 1) * batch_size + 1
        end_idx <- min(i * batch_size, n)
        batch_means[i] <- mean(values[start_idx:end_idx])
      }
      
      # Calculate variance ratio
      var_ratio <- stats::var(batch_means) / (stats::var(values) / batch_size)
      ess <- n / var_ratio
      
      return(list(ess = ess, acf_sum = var_ratio - 1, var_ratio = var_ratio))
      
    } else if (method == "spectral") {
      # Spectral density method (simplified)
      # This is a basic implementation - for production use, consider using coda package
      max_lag <- min(n - 1, floor(n/4))
      acf_result <- stats::acf(values, lag.max = max_lag, plot = FALSE)
      acf_values <- acf_result$acf[, 1, 1]
      
      # Approximate spectral density at frequency 0
      spectral_density <- 1 + 2 * sum(acf_values[-1])
      ess <- n / spectral_density
      
      return(list(ess = ess, acf_sum = spectral_density - 1, var_ratio = spectral_density))
    }
  }
  
  # Extract parameter values based on parameter type
  if (parameter == "all") {
    # Calculate ESS for all available parameters
    available_params <- c("logpost", "kstar")

    # Only pass ktrunc if it is available
    if (all(sapply(samples, function(s) "ktrunc" %in% names(s)))) {
      available_params <- c(available_params, "ktrunc")
    } else {
      warning("ktrunc not available in MCMC samples; skipping it")
    }

    # Check if mu is available
    if (all(sapply(samples, function(s) "mu" %in% names(s)))) {
      available_params <- c(available_params, "mu")
    }
    
    # Check if A and D are available
    if (all(sapply(samples, function(s) "A" %in% names(s)))) {
      available_params <- c(available_params, "A")
    }
    if (all(sapply(samples, function(s) "D" %in% names(s)))) {
      available_params <- c(available_params, "D")
    }
    
    results_list <- list()
    for (param in available_params) {
      results_list[[param]] <- effective_sample_size(results, param, method, additional_burnin)
    }
    
    class(results_list) <- "ess_all_results"
    return(results_list)
    
  } else if (parameter == "logpost") {
    values <- sapply(samples, function(s) s$logpost)
    
  } else if (parameter == "kstar") {
    values <- sapply(samples, function(s) s$kstar)
    
  } else if (parameter == "ktrunc") {
    if (!all(sapply(samples, function(s) "ktrunc" %in% names(s)))) {
      stop("ktrunc parameter not available in MCMC samples")
    }
    values <- sapply(samples, function(s) s$ktrunc)

  } else if (parameter == "n_strains") {
    values <- sapply(samples, function(s) sum(colSums(s$A) > 0))
    
  } else if (parameter == "mu") {
    # Handle stick-breaking weights
    if (!all(sapply(samples, function(s) "mu" %in% names(s)))) {
      stop("mu parameter not available in MCMC samples")
    }
    
    # Get the maximum length of mu across all samples
    max_mu_length <- max(sapply(samples, function(s) length(s$mu)))
    
    # Create matrix of mu values (pad with NA if necessary)
    mu_matrix <- matrix(NA, nrow = n_samples, ncol = max_mu_length)
    for (i in 1:n_samples) {
      mu_length <- length(samples[[i]]$mu)
      mu_matrix[i, 1:mu_length] <- samples[[i]]$mu
    }
    
    # Calculate ESS for each mu component
    ess_results <- list()
    for (j in 1:max_mu_length) {
      mu_values <- mu_matrix[, j]
      # Remove NA values
      valid_values <- mu_values[!is.na(mu_values)]
      if (length(valid_values) > 1) {
        ess_results[[paste0("mu_", j)]] <- calculate_ess(valid_values, method)
        ess_results[[paste0("mu_", j)]]$n_samples <- length(valid_values)
      } else {
        ess_results[[paste0("mu_", j)]] <- list(ess = length(valid_values), 
                                               acf_sum = 0, var_ratio = 1, 
                                               n_samples = length(valid_values))
      }
    }
    
    result <- list(
      parameter = parameter,
      method = method,
      components = ess_results,
      n_samples = n_samples
    )
    
    class(result) <- "ess_result"
    return(result)
    
  } else if (parameter == "A") {
    # Handle allocation matrix
    if (!all(sapply(samples, function(s) "A" %in% names(s)))) {
      stop("A parameter not available in MCMC samples")
    }
    
    # Calculate ESS for each element of A matrix
    A_dim <- dim(samples[[1]]$A)
    ess_results <- list()
    
    for (i in 1:A_dim[1]) {
      for (j in 1:A_dim[2]) {
        A_values <- sapply(samples, function(s) s$A[i, j])
        ess_results[[paste0("A[", i, ",", j, "]")]] <- calculate_ess(A_values, method)
        ess_results[[paste0("A[", i, ",", j, "]")]]$n_samples <- n_samples
      }
    }
    
    # Calculate overall ESS (mean across all elements)
    ess_values <- sapply(ess_results, function(x) x$ess)
    overall_ess <- mean(ess_values, na.rm = TRUE)
    
    result <- list(
      parameter = parameter,
      method = method,
      overall_ess = overall_ess,
      components = ess_results,
      n_samples = n_samples,
      matrix_dimensions = A_dim
    )
    
    class(result) <- "ess_result"
    return(result)
    
  } else if (parameter == "D") {
    # Handle dictionary matrix
    if (!all(sapply(samples, function(s) "D" %in% names(s)))) {
      stop("D parameter not available in MCMC samples")
    }
    
    # Calculate ESS for each element of D matrix
    D_dim <- dim(samples[[1]]$D)
    ess_results <- list()
    
    for (i in 1:D_dim[1]) {
      for (j in 1:D_dim[2]) {
        D_values <- sapply(samples, function(s) s$D[i, j])
        ess_results[[paste0("D[", i, ",", j, "]")]] <- calculate_ess(D_values, method)
        ess_results[[paste0("D[", i, ",", j, "]")]]$n_samples <- n_samples
      }
    }
    
    # Calculate overall ESS (mean across all elements)
    ess_values <- sapply(ess_results, function(x) x$ess)
    overall_ess <- mean(ess_values, na.rm = TRUE)
    
    result <- list(
      parameter = parameter,
      method = method,
      overall_ess = overall_ess,
      components = ess_results,
      n_samples = n_samples,
      matrix_dimensions = D_dim
    )
    
    class(result) <- "ess_result"
    return(result)
    
  } else {
    stop("Invalid parameter. Choose from: 'logpost', 'kstar', 'ktrunc', 'n_strains', 'mu', 'A', 'D', 'all'")
  }
  
  # Calculate ESS for scalar parameters
  ess_result <- calculate_ess(values, method)

  # Add diagnostic information
  result <- list(
    parameter = parameter,
    method = method,
    ess = ess_result$ess,
    n_samples = n_samples,
    acf_sum = ess_result$acf_sum,
    var_ratio = ess_result$var_ratio,
    efficiency = ess_result$ess / n_samples,  # ESS as fraction of total samples
    values = values  # Include actual values for potential further analysis
  )
  
  # Add appropriate class for printing
  class(result) <- "ess_result"
  
  return(result)
}

#' Print summary of SNP-Slice results
#'
#' @param object SNP-Slice results object
#' @param ... Additional arguments
#'
#' @return Summary information
#' @export
#' @importFrom stats median
summary.snp_slice_results <- function(object, ...) {
  if (!inherits(object, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  
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
  cat("- Number of strains identified:", nrow(object$dictionary_matrix), "\n")
  cat("- Number of hosts:", nrow(object$allocation_matrix), "\n")
  
  # Multiplicity of infection
  moi <- rowSums(object$allocation_matrix)
  cat("- Multiplicity of infection (MOI):\n")
  cat("  - Mean MOI:", round(mean(moi), 2), "\n")
  cat("  - Median MOI:", round(stats::median(moi), 2), "\n")
  cat("  - Range:", min(moi), "-", max(moi), "\n")
  
  # Single vs mixed infections
  single_infections <- sum(moi == 1)
  mixed_infections <- sum(moi > 1)
  cat("  - Single infections:", single_infections, "(", round(100 * single_infections / length(moi), 1), "%)\n")
  cat("  - Mixed infections:", mixed_infections, "(", round(100 * mixed_infections / length(moi), 1), "%)\n\n")
  
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
#' @param ... Additional arguments
#'
#' @return Print information
#' @export
print.snp_slice_results <- function(x, ...) {
  if (!inherits(x, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  
  cat("SNP-Slice Results\n")
  cat("================\n")
  cat("Model:", x$model_info$model, "\n")
  cat("Dimensions:", nrow(x$allocation_matrix), "hosts x", ncol(x$dictionary_matrix), "strains x", ncol(x$dictionary_matrix), "SNPs\n")
  cat("Strains identified:", nrow(x$dictionary_matrix), "\n")
  
  if (!is.null(x$convergence)) {
    cat("Gap Converged:", ifelse(x$convergence$gap_converged, "Yes", "No"), "\n")
  }
  
  cat("\n")
}

#' Print effective sample size results
#'
#' @param x ESS result object
#' @param digits Number of digits to display
#' @param ... Additional arguments
#'
#' @return Formatted ESS results
#' @export
print.ess_result <- function(x, digits = 2, ...) {
  cat("Effective Sample Size Results\n")
  cat("============================\n")
  cat("Parameter:", x$parameter, "\n")
  cat("Method:", x$method, "\n")
  cat("Total samples:", x$n_samples, "\n")
  
  if (x$parameter %in% c("logpost", "kstar", "ktrunc", "n_strains")) {
    cat("Effective sample size:", round(x$ess, digits), "\n")
    cat("Efficiency:", round(x$efficiency * 100, 1), "%\n")
    cat("Autocorrelation sum:", round(x$acf_sum, 3), "\n")
    cat("Variance ratio:", round(x$var_ratio, 3), "\n")
  } else if (x$parameter == "mu") {
    cat("ESS for each stick-breaking weight:\n")
    for (name in names(x$components)) {
      comp <- x$components[[name]]
      cat("  ", name, ": ", round(comp$ess, digits), 
          " (", round(comp$ess/comp$n_samples * 100, 1), "% efficiency)\n", sep = "")
    }
  } else if (x$parameter %in% c("A", "D")) {
    cat("Matrix dimensions:", paste(x$matrix_dimensions, collapse = " x "), "\n")
    cat("Overall ESS:", round(x$overall_ess, digits), "\n")
    cat("Overall efficiency:", round(x$overall_ess/x$n_samples * 100, 1), "%\n")
    
    # Show summary statistics for individual elements
    ess_values <- sapply(x$components, function(comp) comp$ess)
    cat("ESS summary:\n")
    cat("  Min:", round(min(ess_values), digits), "\n")
    cat("  Median:", round(median(ess_values), digits), "\n")
    cat("  Max:", round(max(ess_values), digits), "\n")
  }
  cat("\n")
}

#' Print comprehensive ESS results for all parameters
#'
#' @param x List of ESS results
#' @param digits Number of digits to display
#' @param ... Additional arguments
#'
#' @return Formatted ESS results
#' @export
print.ess_all_results <- function(x, digits = 2, ...) {
  cat("Comprehensive Effective Sample Size Analysis\n")
  cat("===========================================\n\n")
  
  for (param_name in names(x)) {
    cat("Parameter:", param_name, "\n")
    cat(paste(rep("-", 12 + nchar(param_name)), collapse = ""), "\n")
    print(x[[param_name]], digits = digits)
  }
}
