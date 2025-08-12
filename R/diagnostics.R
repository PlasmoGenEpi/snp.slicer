# Global variables for data.frame column names
utils::globalVariables(c("iteration", "logpost", "kstar", "n_strains"))

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
  allocations <- list(
    allocation_matrix = results$allocation_matrix,
    n_hosts = nrow(results$allocation_matrix),
    n_strains = ncol(results$allocation_matrix),
    host_names = paste0("Host_", 1:nrow(results$allocation_matrix)),
    strain_names = paste0("Strain_", 1:ncol(results$allocation_matrix)),
    multiplicity_of_infection = rowSums(results$allocation_matrix)
  )
  
  return(allocations)
}

#' Plot convergence diagnostics
#'
#' @param results SNP-Slice results object
#' @param type Type of plot ("logpost", "kstar", "n_strains")
#'
#' @return Plot object (if ggplot2 is available)
#' @export
plot_convergence <- function(results, type = "logpost") {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  
  if (is.null(results$mcmc_samples)) {
    stop("MCMC samples not stored. Set store_mcmc = TRUE when running snp_slice()")
  }
  
  # Extract MCMC samples
  samples <- results$mcmc_samples
  
  if (type == "logpost") {
    logpost_values <- sapply(samples, function(s) s$logpost)
    iterations <- 1:length(logpost_values)
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
    iterations <- 1:length(kstar_values)
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
    iterations <- 1:length(n_strains_values)
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
  } else {
    stop("Invalid plot type. Choose from: 'logpost', 'kstar', 'n_strains'")
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
#'
#' @return Effective sample size(s) and diagnostic information
#' @export
#' @importFrom stats acf var
effective_sample_size <- function(results, parameter = "logpost", method = "autocorrelation") {
  if (!inherits(results, "snp_slice_results")) {
    stop("Input must be an snp_slice_results object")
  }
  
  if (is.null(results$mcmc_samples)) {
    stop("MCMC samples not stored. Set store_mcmc = TRUE when running snp_slice()")
  }
  
  samples <- results$mcmc_samples
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
    available_params <- c("logpost", "kstar", "ktrunc")
    
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
      results_list[[param]] <- effective_sample_size(results, param, method)
    }
    
    class(results_list) <- "ess_all_results"
    return(results_list)
    
  } else if (parameter == "logpost") {
    values <- sapply(samples, function(s) s$logpost)
    
  } else if (parameter == "kstar") {
    values <- sapply(samples, function(s) s$kstar)
    
  } else if (parameter == "ktrunc") {
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
  cat("Model:", object$model_info$model, "\n")
  cat("Data dimensions:", object$model_info$N, "hosts x", object$model_info$P, "SNPs\n")
  cat("Data type:", object$model_info$data_type, "\n\n")
  
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
    cat("- Converged:", ifelse(object$convergence$converged, "Yes", "No"), "\n")
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
    cat("Converged:", ifelse(x$convergence$converged, "Yes", "No"), "\n")
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
