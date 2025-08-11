#' Create model object
#'
#' @param model Model type
#' @param processed_data Processed data
#' @param ... Additional model parameters
#'
#' @return Model object with initialization and likelihood functions
#' @keywords internal
create_model <- function(model, processed_data, alpha = 2.6, rho = 0.5, ...) {
  
  # Extract additional parameters
  params <- list(...)
  
  # Create model-specific object
  if (model == "categorical") {
    e1 <- or_null(params$e1, 0.05)
    e2 <- or_null(params$e2, 0.05)
    
    model_obj <- list(
      name = "categorical",
      y = processed_data$y,
      r = processed_data$r,
      N = processed_data$N,
      P = processed_data$P,
      alpha = alpha,
      rho = rho,
      e1 = e1,
      e2 = e2,
      loglikelihood_matrix = categorical_loglikelihood_matrix,
      loglikelihood_vector = categorical_loglikelihood_vector,
      initialize_state = categorical_initialize_state,
      resolve_exceptions = categorical_resolve_exceptions
    )
  } else if (model == "poisson") {
    model_obj <- list(
      name = "poisson",
      y = processed_data$y,
      r = processed_data$r,
      N = processed_data$N,
      P = processed_data$P,
      alpha = alpha,
      rho = rho,
      loglikelihood_matrix = poisson_loglikelihood_matrix,
      loglikelihood_vector = poisson_loglikelihood_vector,
      initialize_state = poisson_initialize_state,
      resolve_exceptions = poisson_resolve_exceptions
    )
  } else if (model == "binomial") {
    model_obj <- list(
      name = "binomial",
      y = processed_data$y,
      r = processed_data$r,
      N = processed_data$N,
      P = processed_data$P,
      alpha = alpha,
      rho = rho,
      loglikelihood_matrix = binomial_loglikelihood_matrix,
      loglikelihood_vector = binomial_loglikelihood_vector,
      initialize_state = binomial_initialize_state,
      resolve_exceptions = binomial_resolve_exceptions
    )
  } else if (model == "negative_binomial") {
    model_obj <- list(
      name = "negative_binomial",
      y = processed_data$y,
      r = processed_data$r,
      N = processed_data$N,
      P = processed_data$P,
      alpha = alpha,
      rho = rho,
      loglikelihood_matrix = negative_binomial_loglikelihood_matrix,
      loglikelihood_vector = negative_binomial_loglikelihood_vector,
      initialize_state = negative_binomial_initialize_state,
      resolve_exceptions = negative_binomial_resolve_exceptions
    )
  } else {
    stop("Unsupported model: ", model)
  }
  
  class(model_obj) <- "snp_slice_model"
  return(model_obj)
}

#' Run MCMC for SNP-Slice
#'
#' @param model_obj Model object
#' @param n_mcmc Number of MCMC iterations
#' @param burnin Burn-in period
#' @param gap Early stopping threshold
#' @param verbose Whether to print progress
#' @param store_mcmc Whether to store full MCMC samples
#'
#' @return MCMC results
#' @keywords internal
run_mcmc <- function(model_obj, n_mcmc, burnin, gap, verbose, store_mcmc) {
  
  # Initialize state
  if (verbose) {
    cat("Initializing state...\n")
  }
  
  state <- model_obj$initialize_state(model_obj, threshold = 0.001)
  
  if (verbose) {
    cat("Initialization complete, log-likelihood:", state$loglik, "\n")
    cat("Starting with", sum(rowSums(state$A) == 1), "single infections\n")
  }
  
  # Initialize slice sampler
  state <- slice_init(state, model_obj)
  
  if (verbose) {
    cat("Plan to run", n_mcmc, "iterations with", burnin, "burnin, gap =", gap, "\n")
  }
  
  # Clear performance log before starting
  clear_performance_log()
  
  # Run MCMC
  map_state <- state
  lpostmax <- -Inf
  mapiter <- 0
  mapktrunc <- state$ktrunc
  
  # Storage for samples if requested
  samples <- if (store_mcmc) list() else NULL
  
  for (iter in 1:n_mcmc) {
    # Run one MCMC iteration
    state <- slice_iter(state, model_obj)
    
    # Store sample if requested
    if (store_mcmc) {
      samples[[iter]] <- list(
        A = state$A,
        D = state$D,
        mu = state$mu,
        logpost = state$logpost,
        kstar = state$kstar
      )
    }
    
    # Update MAP estimate
    if (state$ktrunc > mapktrunc) {
      # If ktrunc changes, restart MAP
      map_state <- state
      mapiter <- iter
      mapktrunc <- state$ktrunc
      lpostmax <- state$logpost
    } else if (state$logpost > lpostmax) {
      # If same ktrunc but higher logpost
      map_state <- state
      mapiter <- iter
      lpostmax <- state$logpost
    }
    
    # Print progress
    if (verbose && iter %% 100 == 0) {
      cat("Iteration", iter, "active strains:", sum(colSums(state$A) > 0),
          "kstar:", state$kstar, "ktrunc:", state$ktrunc,
          "single infections:", sum(rowSums(state$A) == 1),
          "logpost:", round(state$logpost, 2), "max:", round(lpostmax, 2), "\n")
    }
    
    # Check for early stopping
    if (iter > burnin && !is.null(gap) && mapiter < iter - gap) {
      if (verbose) {
        cat("Early stopping at iteration", iter, "due to convergence\n")
      }
      break
    }
  }
  
  # Create diagnostics
  diagnostics <- list(
    final_iteration = iter,
    map_iteration = mapiter,
    final_logpost = state$logpost,
    map_logpost = lpostmax,
    final_kstar = state$kstar,
    map_kstar = map_state$kstar,
    final_ktrunc = state$ktrunc,
    map_ktrunc = map_state$ktrunc
  )
  
  # Create results
  result <- list(
    map_state = map_state,
    final_state = state,
    samples = samples,
    diagnostics = diagnostics,
    parameters = list(
      n_mcmc = n_mcmc,
      burnin = burnin,
      gap = gap,
      store_mcmc = store_mcmc
    ),
    convergence = list(
      converged = !is.null(gap) && mapiter < iter - gap,
      iterations_run = iter
    )
  )
  
  return(result)
}

#' Initialize slice sampler
#'
#' @param state Initial state
#' @param model_obj Model object
#'
#' @return State with slice sampler initialized
#' @keywords internal
#' @importFrom stats rbeta runif
slice_init <- function(state, model_obj) {
  
  # Initialize stick-breaking weights
  ktrunc <- ncol(state$A)
  mu <- stats::rbeta(ktrunc + 1, model_obj$alpha / (ktrunc + 1), 1)
  mu_sort <- sort(mu, decreasing = TRUE)
  
  # Update state
  state$mu <- mu_sort
  state$ktrunc <- ktrunc
  state$kstar <- get_kstar(state$A)
  state$kplus <- state$kstar + 1
  
  # Expand matrices
  state$A <- cbind(state$A, 0)
  state$D <- rbind(state$D, stats::runif(model_obj$P) < model_obj$rho)
  state$ktrunc <- ncol(state$A)
  
  return(state)
}

#' Get last active feature index
#'
#' @param A Allocation matrix
#'
#' @return Index of last active feature
#' @keywords internal
get_kstar <- function(A) {
  active_cols <- which(colSums(A) > 0)
  if (length(active_cols) == 0) {
    return(0)
  }
  return(max(active_cols))
}

#' Get mu value for last active feature
#'
#' @param A Allocation matrix
#' @param mu Stick-breaking weights
#'
#' @return Mu value for last active feature
#' @keywords internal
get_mustar <- function(A, mu) {
  kstar <- get_kstar(A)
  if (kstar == 0) {
    return(0)
  }
  return(mu[kstar])
}

#' Get column sums of allocation matrix
#'
#' @param A Allocation matrix
#'
#' @return Column sums
#' @keywords internal
get_m <- function(A) {
  return(colSums(A))
}

#' Remove duplicates from strain patterns
#'
#' @param dd Binary matrix of strain patterns
#'
#' @return List with assignments and unique dictionary
#' @keywords internal
remove_duplicates <- function(dd) {
  # Calculate pairwise distances
  metric <- 1 - (dd %*% t(dd) + (1 - dd) %*% t(1 - dd)) / ncol(dd)
  
  # Find unique patterns
  assignments <- rep(NA, nrow(dd))
  counter <- 0
  first_appear <- c()
  
  for (i in 1:nrow(dd)) {
    if (is.na(assignments[i])) {
      counter <- counter + 1
      twins <- which(metric[i, ] == 0)
      assignments[twins] <- counter
      first_appear <- c(first_appear, i)
    }
  }
  
  return(list(
    assignments = assignments,
    D = dd[first_appear, , drop = FALSE]
  ))
}

#' Helper function for default values
#'
#' @param x Value to check
#' @param default Default value
#'
#' @return x if not NULL, otherwise default
#' @keywords internal
or_null <- function(x, default) {
  if (is.null(x)) default else x
}

#' Null coalescing operator (infix)
#'
#' @param x Value to check
#' @param default Default value if x is NULL
#'
#' @return x if not NULL, otherwise default
#' @keywords internal
`%||%` <- function(x, default) {
  if (is.null(x)) default else x
}
