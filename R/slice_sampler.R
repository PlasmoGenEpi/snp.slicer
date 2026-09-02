#' Run one slice sampling iteration
#'
#' @param state Current state
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
#' @importFrom stats runif
slice_iter <- function(state, model_obj) {
  start_timer("slice_iter")

  kernel <- model_obj$kernel
  if (is.null(kernel)) {
    kernel <- mcmc_kernel_r()
  }

  if (identical(kernel$name, "cpp") && is.null(kernel$obs_code)) {
    kernel <- kernel_with_obs_cache(kernel, model_obj)
    model_obj$kernel <- kernel
    model_obj$kernel_loglik_const <- kernel$obs_loglik_const
    model_obj$kernel_obs_code <- kernel$obs_code
  }

  if (!is.null(kernel$update_iter)) {
    start_timer("slice_update_iter")
    state <- kernel$update_iter(state, model_obj)
    end_timer("slice_update_iter")
  } else {
    state <- slice_update_s(state, model_obj)
    state <- slice_update_a(state, model_obj)
    state <- slice_update_d(state, model_obj)
    state <- slice_update_mu(state, model_obj)
  }

  # Update log-likelihood and log-posterior
  if (identical(kernel$name, "cpp") && exists("cpp_loglik_total", where = asNamespace("snp.slicer"), inherits = FALSE)) {
    obs <- kernel_obs_args(model_obj)
    state$loglik <- as.numeric(cpp_loglik_total(
      A = state$A,
      D = state$D,
      y = model_obj$y,
      r = model_obj$r,
      model_type = model_type_id(model_obj$name),
      llik_tab = kernel_llik_tab(model_obj),
      loglik_const = obs$loglik_const,
      obs_code = obs$obs_code
    ))
  } else {
    state$loglik <- model_obj$loglikelihood_matrix(state$A, state$D, model_obj)
  }
  state$logpost <- as.numeric(state$loglik) +
                   logprior_a(state$A, state$mu, model_obj$alpha) +
                   logprior_mu(state$mu, model_obj$alpha, model_obj$N) +
                   logprior_d(state$D, model_obj$rho)
  
  end_timer("slice_iter")
  return(state)
}

#' Update slice variable s
#'
#' @param state Current state
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
#' @importFrom stats runif
slice_update_s <- function(state, model_obj) {
  start_timer("slice_update_s")
  kernel <- model_obj$kernel
  if (is.null(kernel)) {
    kernel <- mcmc_kernel_r()
  }
  state <- kernel$update_s(state, model_obj)
  end_timer("slice_update_s")
  return(state)
}

#' R reference implementation of slice-variable update
#' @keywords internal
slice_update_s_r <- function(state, model_obj) {
  mustar <- get_mustar(state$A, state$mu)
  s <- mustar * stats::runif(1)

  k <- state$ktrunc
  while (s < state$mu[k]) {
    munext <- gridsample(0, state$mu[k], logf_newfeature, model_obj$N, model_obj$alpha)
    state$mu <- c(state$mu, munext)
    k <- k + 1
  }

  if (state$ktrunc < length(state$mu)) {
    num_new_cols <- length(state$mu) - state$ktrunc
    state$A <- cbind(state$A, matrix(0, nrow = model_obj$N, ncol = num_new_cols))
    state$D <- rbind(state$D, matrix(0, ncol = model_obj$P, nrow = num_new_cols))

    for (k in (state$ktrunc + 1):length(state$mu)) {
      state <- refresh_feature(state, k, model_obj$rho, model_obj$P)
    }
  }

  state$ktrunc <- length(state$mu)
  state$kplus <- which(state$mu < s)[1]
  if (is.na(state$kplus)) {
    state$kplus <- state$ktrunc
  }

  return(state)
}

#' Update allocation matrix A
#'
#' @param state Current state
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
slice_update_a <- function(state, model_obj) {
  start_timer("slice_update_a")
  kernel <- model_obj$kernel
  if (is.null(kernel)) {
    kernel <- mcmc_kernel_r()
  }
  state <- kernel$update_a(state, model_obj)
  end_timer("slice_update_a")
  return(state)
}

#' R reference implementation of allocation-matrix updates
#' @keywords internal
slice_update_a_r <- function(state, model_obj) {
  kstar <- state$kstar

  for (i in state$mixed) {
    full_row_i <- state$A[i, , drop = FALSE] %*% state$D
    row_sum_i <- sum(state$A[i, ])
    for (k in 1:state$kplus) {
      old_A_ik <- state$A[i, k]
      state <- slice_update_a_local(state, i, k, model_obj,
        full_row_i = full_row_i, row_sum_i = row_sum_i, kstar = kstar)
      if (state$A[i, k] != old_A_ik) {
        delta <- state$A[i, k] - old_A_ik
        full_row_i <- full_row_i + delta * state$D[k, ]
        row_sum_i <- row_sum_i + delta
        if (state$A[i, k] == 1L) {
          kstar <- max(kstar, k)
        } else if (k == kstar && sum(state$A[, k]) == 0) {
          active <- which(colSums(state$A) > 0)
          kstar <- if (length(active) == 0) 0 else max(active)
        }
      }
    }
  }
  state$kstar <- kstar
  return(state)
}

#' Update single element of allocation matrix A
#'
#' @param state Current state
#' @param i Host index
#' @param k Strain index
#' @param model_obj Model object
#' @param full_row_i Optional precomputed row contribution \code{A[i,]} %*% D (from slice_update_a)
#' @param row_sum_i Optional precomputed sum of \code{A[i,]} (from slice_update_a)
#' @param kstar Optional current last-active-feature index (from slice_update_a)
#'
#' @return Updated state
#' @keywords internal
#' @importFrom stats runif
slice_update_a_local <- function(state, i, k, model_obj,
                                 full_row_i = NULL, row_sum_i = NULL, kstar = NULL) {
  start_timer("slice_update_a_local")

  if (is.null(full_row_i) || is.null(row_sum_i)) {
    ad0 <- state$A[i, -k, drop = FALSE] %*% state$D[-k, , drop = FALSE]
    a0 <- sum(state$A[i, -k])
  } else {
    ad0 <- full_row_i - state$A[i, k] * state$D[k, ]
    a0 <- row_sum_i - state$A[i, k]
  }
  if (is.null(kstar)) {
    kstar <- state$kstar
  }

  if (a0 == 0) {
    state$A[i, k] <- 1
    return(state)
  }

  # Calculate log probabilities
  logp0 <- model_obj$loglikelihood_vector(as.vector(ad0 / a0),
                                         as.vector(model_obj$y[i, ]),
                                         as.vector(model_obj$r[i, ]))
  logp1 <- model_obj$loglikelihood_vector(as.vector((ad0 + state$D[k, ]) / (a0 + 1)),
                                         as.vector(model_obj$y[i, ]),
                                         as.vector(model_obj$r[i, ]))

  # Add prior contributions
  logp0 <- logp0 + log(1 - state$mu[k])
  logp1 <- logp1 + log(state$mu[k])

  # Handle special cases for kstar changes (use passed-in kstar)
  if (k == kstar && state$A[i, k] == 1 && sum(state$A[-i, k]) == 0) {
    logp1 <- logp1 - log(state$mu[k])
    next_kstar <- tail(which(colSums(state$A) != 0), 2)[1]
    logp0 <- logp0 - log(state$mu[next_kstar])
  } else if (k > kstar) {
    logp1 <- logp1 - log(state$mu[k])
    logp0 <- logp0 - log(state$mu[kstar])
  }

  # Metropolis-Hastings acceptance
  p1 <- get_mhratio(logp1, logp0)
  u <- stats::runif(1)
  state$A[i, k] <- (u < p1)

  end_timer("slice_update_a_local")
  return(state)
}

#' Update dictionary matrix D
#'
#' @param state Current state
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
slice_update_d <- function(state, model_obj) {
  start_timer("slice_update_d")
  kernel <- model_obj$kernel
  if (is.null(kernel)) {
    kernel <- mcmc_kernel_r()
  }
  state <- kernel$update_d(state, model_obj)
  end_timer("slice_update_d")
  return(state)
}

#' R reference implementation of dictionary-matrix updates
#' @keywords internal
slice_update_d_r <- function(state, model_obj) {
  an <- rowSums(state$A)

  for (k in state$kmin:state$kstar) {
    # Leave-k-out product: one matrix multiply per k instead of P per k
    M_k <- state$A[, -k, drop = FALSE] %*% state$D[-k, , drop = FALSE]
    for (p in 1:model_obj$P) {
      state <- slice_update_d_local(state, k, p, model_obj, an = an, ad0 = M_k[, p])
    }
  }

  return(state)
}

#' Update single element of dictionary matrix D
#'
#' @param state Current state
#' @param k Strain index
#' @param p SNP index
#' @param model_obj Model object
#' @param an Optional precomputed row sums of A (from slice_update_d)
#' @param ad0 Optional precomputed leave-k-out contribution for this (k,p)
#'
#' @return Updated state
#' @keywords internal
#' @importFrom stats runif
slice_update_d_local <- function(state, k, p, model_obj, an = NULL, ad0 = NULL) {
  start_timer("slice_update_d_local")

  if (is.null(an)) {
    an <- rowSums(state$A)
  }
  if (is.null(ad0)) {
    ad0 <- state$A[, -k, drop = FALSE] %*% state$D[-k, p, drop = FALSE]
  }

  # Prior probabilities
  logp1 <- log(model_obj$rho)
  logp0 <- log(1 - model_obj$rho)

  # Likelihood contributions
  logp1 <- logp1 + model_obj$loglikelihood_vector((ad0 + state$A[, k]) / an,
                                                 model_obj$y[, p],
                                                 model_obj$r[, p])
  logp0 <- logp0 + model_obj$loglikelihood_vector(ad0 / an,
                                                 model_obj$y[, p],
                                                 model_obj$r[, p])

  # Metropolis-Hastings acceptance
  p1 <- get_mhratio(logp1, logp0)
  u <- stats::runif(1)
  state$D[k, p] <- (u < p1)

  end_timer("slice_update_d_local")
  return(state)
}

#' Update stick-breaking weights mu
#'
#' @param state Current state
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
slice_update_mu <- function(state, model_obj) {
  start_timer("slice_update_mu")
  kernel <- model_obj$kernel
  if (is.null(kernel)) {
    kernel <- mcmc_kernel_r()
  }
  state <- kernel$update_mu(state, model_obj)
  end_timer("slice_update_mu")
  return(state)
}

#' R reference implementation of stick-breaking weight updates
#' @keywords internal
slice_update_mu_r <- function(state, model_obj) {
  m <- get_m(state$A)
  
  # Update mu[1] (first stick)
  if (state$kplus > 1) {
    if (2 <= length(state$mu) && !is.na(state$mu[2]) && is.finite(state$mu[2])) {
      state$mu[1] <- gridsample(state$mu[2], 1, function(x) logf_oldfeature(x, m[1], model_obj$N))
    } else {
      state$mu[1] <- 0.5  # Default value
    }
  }
  
  # Update mu[2:kplus-1] (middle sticks)
  if (state$kplus > 2L) {
    for (k in 2:(state$kplus - 1)) {
      if (k + 1 <= length(state$mu) && k - 1 >= 1 && 
          !is.na(state$mu[k + 1]) && !is.na(state$mu[k - 1]) && 
          is.finite(state$mu[k + 1]) && is.finite(state$mu[k - 1])) {
        state$mu[k] <- gridsample(state$mu[k + 1], state$mu[k - 1], function(x) logf_oldfeature(x, m[k], model_obj$N))
      } else {
        state$mu[k] <- 0.5  # Default value
      }
    }
  }
  
  # Update mu[kplus:ktrunc] (inactive features)
  for (k in state$kplus:state$ktrunc) {
    if (k - 1 >= 1 && !is.na(state$mu[k - 1]) && is.finite(state$mu[k - 1])) {
      state$mu[k] <- gridsample(0, state$mu[k - 1], logf_newfeature, model_obj$N, model_obj$alpha)
    } else {
      state$mu[k] <- 0.5  # Default value
    }
  }

  return(state)
}

#' Metropolis-Hastings acceptance ratio
#'
#' @param logp1 Log probability of state 1
#' @param logp0 Log probability of state 0
#'
#' @return Acceptance probability
#' @keywords internal
get_mhratio <- function(logp1, logp0) {
  if (is.infinite(logp1)) return(0)
  if (is.infinite(logp0)) return(1)
  
  maxlogp <- max(logp1, logp0)
  mhratio <- exp(logp1 - maxlogp) / (exp(logp1 - maxlogp) + exp(logp0 - maxlogp))
  return(mhratio)
}

#' Grid sampling for constrained updates
#'
#' @param lb Lower bound
#' @param ub Upper bound
#' @param logf Log density function
#' @param ... Additional arguments for logf
#'
#' @return Sampled value
#' @keywords internal
gridsample <- function(lb, ub, logf, ...) {
  start_timer("gridsample")
  
  # Check for valid bounds
  if (!is.finite(lb) || !is.finite(ub)) {
    stop("Grid sampling requires finite bounds")
  }
  
  if (lb >= ub) {
    # If bounds are invalid, return a reasonable default
    return(0.5)
  }
  
  nsample <- 100
  xgrid <- seq(lb, ub, length.out = nsample)
  xgrid <- xgrid[-c(1, nsample)]  # Exclude endpoints
  
  lw <- logf(xgrid, ...)
  result <- sample_index_logweights(lw, xgrid)
  
  end_timer("gridsample")
  return(result)
}

#' Sample from categorical distribution with log weights
#'
#' @param lw Log weights
#' @param xgrid Grid of values
#'
#' @return Sampled value
#' @keywords internal
sample_index_logweights <- function(lw, xgrid) {
  lw <- lw - max(lw)
  probs <- exp(lw) / sum(exp(lw))
  idx <- sample.int(length(lw), size = 1, prob = probs)
  return(xgrid[idx])
}

#' Log density for new features
#'
#' @param x Value
#' @param N Number of hosts
#' @param alpha IBP concentration parameter
#'
#' @return Log density
#' @keywords internal
logf_newfeature <- function(x, N, alpha) {
  start_timer("logf_newfeature")
  
  logf <- alpha * sum(sapply(1:N, function(i) (1 - x)^i / i))
  logf <- logf + (alpha - 1) * log(x) + N * log(1 - x)
  
  end_timer("logf_newfeature")
  return(logf)
}

#' Log density for existing features
#'
#' @param x Value
#' @param m Feature count
#' @param N Number of hosts
#'
#' @return Log density
#' @keywords internal
logf_oldfeature <- function(x, m, N) {
  return((m - 1) * log(x) + (N - m) * log(1 - x))
}

#' Refresh feature with random dictionary
#'
#' @param state Current state
#' @param k Feature index
#' @param rho Dictionary sparsity parameter
#' @param P Number of SNPs
#'
#' @return Updated state
#' @keywords internal
#' @importFrom stats runif
refresh_feature <- function(state, k, rho, P) {
  state$A[, k] <- 0
  state$D[k, ] <- (stats::runif(P) < rho)
  return(state)
}
