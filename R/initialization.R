#' Create model object
#'
#' @param model Model type
#' @param processed_data Processed data
#' @param ... Additional model parameters
#'
#' @return Model object with initialization and likelihood functions
#' @keywords internal
create_model <- function(
                         model, 
                         processed_data, 
                         alpha = 2.6,
                         rho = if (model == "categorical") 0.5 else NULL,
                         ...) {
  
  # Extract additional parameters
  params <- list(...)
  if (is.null(rho)) {
    polymorphic_obs <- !is.na(processed_data$y) & processed_data$y != processed_data$r
    rho <- sum(processed_data$y[polymorphic_obs], na.rm = TRUE) /
      sum(processed_data$r[polymorphic_obs], na.rm = TRUE)
  }
  
  # Create model-specific object
  if (model == "categorical") {
    e1 <- or_null(params$e1, 0.05)
    e2 <- or_null(params$e2, 0.05)
    llik_tab <- build_categorical_llik_tab(e1, e2)
    r_mat <- processed_data$r
    if (is.null(r_mat)) {
      r_mat <- matrix(0, nrow = processed_data$N, ncol = processed_data$P)
    }

    model_obj <- list(
      name = "categorical",
      y = processed_data$y,
      r = r_mat,
      N = processed_data$N,
      P = processed_data$P,
      alpha = alpha,
      rho = rho,
      e1 = e1,
      e2 = e2,
      llik_tab = llik_tab,
      loglikelihood_matrix = categorical_loglikelihood_matrix,
      loglikelihood_vector = function(propvec, yvec, rvec = NULL) {
        categorical_loglikelihood_vector(propvec, yvec, llik_tab = llik_tab)
      },
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

#' Run a single MCMC chain for SNP-Slice
#'
#' This function runs a single MCMC chain for a configured snp.slicer 
#' model object.
#'
#' @param model_obj Model object
#' @param n_sample Number of post-burn-in iterations to retain
#' @param n_burnin Number of iterations discarded before sampling begins
#' @param gap Early stopping threshold
#' @param verbose Whether to print progress
#' @param store_mcmc Whether to store full MCMC samples
#' @param chain_id Integer identifier of this chain (1 when running one chain)
#' @param seed Integer seed for this chain
#'
#' @details Sampling iterations are returned when \code{store_mcmc = TRUE}.
#'
#' @return MCMC results
#' @keywords internal
run_chain <- function(model_obj, n_sample, n_burnin, gap, verbose, store_mcmc,
                      chain_id = 1L, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Prefix progress output so interleaved chain output stays readable
  prefix <- if (is.null(chain_id)) "" else paste0("[chain ", chain_id, "] ")
  cat_w_prefix <- function(...) cat(prefix, ..., sep = "")

  # Initialize state
  if (verbose) {
    cat_w_prefix("Initializing state...\n")
  }

  state <- model_obj$initialize_state(model_obj, threshold = 0.001)
  
  if (verbose) {
    cat_w_prefix("Initialization complete, log-likelihood: ", state$loglik, "\n")
    cat_w_prefix("Starting with ", sum(rowSums(state$A) == 1), " single infections\n")
  }
  
  # Initialize slice sampler
  state <- slice_init(state, model_obj)
  
  total_iter <- n_burnin + n_sample

  if (verbose) {
    cat_w_prefix("Plan to run ", total_iter, " total iterations: ", n_burnin, " burn-in + ",
        n_sample, " sampling, gap = ", if (is.null(gap)) "NULL" else gap, "\n")
  }

  model_obj <- setup_mcmc_kernel(model_obj)

  # Clear performance log before starting
  clear_performance_log()

  # Run MCMC
  map_state <- state
  lpostmax <- -Inf
  mapiter <- 0
  mapktrunc <- state$ktrunc

  # Storage for post-burn-in samples if requested
  samples <- if (store_mcmc) vector("list", n_sample) else NULL
  n_stored <- 0
  for (iter in 1:total_iter) {
    # Run one MCMC iteration
    state <- slice_iter(state, model_obj)

    # Burn-in iterations are skipped
    post_burnin <- iter > n_burnin

    # Store sample if requested
    if (store_mcmc && post_burnin) {
      n_stored <- n_stored + 1
      samples[[n_stored]] <- list(
        A = state$A,
        D = state$D,
        mu = state$mu,
        logpost = state$logpost,
        kstar = state$kstar,
        ktrunc = state$ktrunc,
        iteration = iter
      )
    }

    # Update MAP estimate
    if (post_burnin) {
      if (mapiter == 0) {
        # Seed the MAP from the first post-burn-in state
        map_state <- state
        mapiter <- iter
        mapktrunc <- state$ktrunc
        lpostmax <- state$logpost
      } else if (state$ktrunc > mapktrunc) {
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
    }

    # Print progress
    if (verbose && iter %% 10 == 0) {
      cat_w_prefix("Iteration ", iter, " active strains: ", sum(colSums(state$A) > 0),
          " kstar: ", state$kstar, " ktrunc: ", state$ktrunc,
          " single infections: ", sum(rowSums(state$A) == 1),
          " logpost: ", round(state$logpost, 2), " max: ", round(lpostmax, 2), "\n")
    }

    # Check for early stopping
    if (iter > n_burnin && !is.null(gap) && mapiter < iter - gap) {
      if (verbose) {
        cat_w_prefix("Early stopping at iteration ", iter, " due to convergence\n")
      }
      break
    }
  }

  # Drop unused slots if the chain stopped early
  if (store_mcmc && n_stored < n_sample) {
    length(samples) <- n_stored
  }

  # Create diagnostics
  diagnostics <- list(
    chain_id = chain_id,
    seed = seed,
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
    chain_id = chain_id,
    seed = seed,
    map_state = map_state,
    final_state = state,
    samples = samples,
    diagnostics = diagnostics,
    parameters = list(
      n_sample = n_sample,
      n_burnin = n_burnin,
      gap = gap,
      store_mcmc = store_mcmc
    ),
    convergence = list(
      gap_converged = !is.null(gap) && mapiter < iter - gap,
      iterations_run = iter,
      samples_retained = max(0, iter - n_burnin)
    ),
    performance = get_performance_summary()
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
  rownames(dd) <- NULL
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


