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
  
  state <- slice_update_s(state, model_obj)
  state <- slice_update_a(state, model_obj)
  state <- slice_update_d(state, model_obj)
  state <- slice_update_mu(state, model_obj)
  
  # Update log-likelihood and log-posterior
  state$loglik <- model_obj$loglikelihood_matrix(state$A, state$D, model_obj)
  state$logpost <- state$loglik + 
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
  
  # Sample auxiliary variable s
  mustar <- get_mustar(state$A, state$mu)
  s <- mustar * stats::runif(1)  # 0 < s < mu[kstar]
  
  # Expand truncation level if needed
  k <- state$ktrunc
  while (s < state$mu[k]) {
    # Sample new mu value using grid sampling
    munext <- gridsample(0, state$mu[k], logf_newfeature, model_obj$N, model_obj$alpha)
    state$mu <- c(state$mu, munext)
    k <- k + 1
  }
  
  # Expand A and D matrices if necessary
  if (state$ktrunc < length(state$mu)) {
    num_new_cols <- length(state$mu) - state$ktrunc
    state$A <- cbind(state$A, matrix(0, nrow = model_obj$N, ncol = num_new_cols))
    state$D <- rbind(state$D, matrix(0, ncol = model_obj$P, nrow = num_new_cols))
    
    # Refresh new features
    for (k in (state$ktrunc + 1):length(state$mu)) {
      state <- refresh_feature(state, k, model_obj$rho, model_obj$P)
    }
  }
  
  state$ktrunc <- length(state$mu)
  state$kplus <- which(state$mu < s)[1]
  if (is.na(state$kplus)) {
    state$kplus <- state$ktrunc
  }
  
  end_timer("slice_update_s")
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
  
  for (k in 1:state$kplus) {
    for (i in state$mixed) {
      state <- slice_update_a_local(state, i, k, model_obj)
    }
  }
  state$kstar <- get_kstar(state$A)
  
  end_timer("slice_update_a")
  return(state)
}

#' Update single element of allocation matrix A
#'
#' @param state Current state
#' @param i Host index
#' @param k Strain index
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
#' @importFrom stats runif
slice_update_a_local <- function(state, i, k, model_obj) {
  start_timer("slice_update_a_local")
  
  # Calculate current allocation excluding strain k
  ad0 <- state$A[i, -k, drop = FALSE] %*% state$D[-k, , drop = FALSE]
  a0 <- sum(state$A[i, -k])
  
  if (a0 == 0) {
    state$A[i, k] <- 1
    state$kstar <- get_kstar(state$A)
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
  
  # Handle special cases for kstar changes
  if (k == state$kstar && state$A[i, k] == 1 && sum(state$A[-i, k]) == 0) {
    logp1 <- logp1 - log(state$mu[k])
    next_kstar <- tail(which(colSums(state$A) != 0), 2)[1]
    logp0 <- logp0 - log(state$mu[next_kstar])
  } else if (k > state$kstar) {
    logp1 <- logp1 - log(state$mu[k])
    logp0 <- logp0 - log(state$mu[state$kstar])
  }
  
  # Metropolis-Hastings acceptance
  p1 <- get_mhratio(logp1, logp0)
  u <- stats::runif(1)
  state$A[i, k] <- (u < p1)
  state$kstar <- get_kstar(state$A)
  
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
  
  for (k in state$kmin:state$kstar) {
    for (p in 1:model_obj$P) {
      state <- slice_update_d_local(state, k, p, model_obj)
    }
  }
  
  end_timer("slice_update_d")
  return(state)
}

#' Update single element of dictionary matrix D
#'
#' @param state Current state
#' @param k Strain index
#' @param p SNP index
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
#' @importFrom stats runif
slice_update_d_local <- function(state, k, p, model_obj) {
  start_timer("slice_update_d_local")
  
  # Calculate current dictionary contribution excluding strain k
  ad0 <- state$A[, -k, drop = FALSE] %*% state$D[-k, p, drop = FALSE]
  an <- rowSums(state$A)
  
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
  for (k in 2:(state$kplus - 1)) {
    if (k + 1 <= length(state$mu) && k - 1 >= 1 && 
        !is.na(state$mu[k + 1]) && !is.na(state$mu[k - 1]) && 
        is.finite(state$mu[k + 1]) && is.finite(state$mu[k - 1])) {
      state$mu[k] <- gridsample(state$mu[k + 1], state$mu[k - 1], function(x) logf_oldfeature(x, m[k], model_obj$N))
    } else {
      state$mu[k] <- 0.5  # Default value
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
  
  end_timer("slice_update_mu")
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
