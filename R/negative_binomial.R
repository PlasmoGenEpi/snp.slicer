#' Negative Binomial Model for SNP-Slice
#'
#' @description
#' Implementation of the negative binomial observation model for SNP-Slice.
#' This is the recommended default model for most applications.

#' Log-likelihood for negative binomial model (matrix version)
#'
#' @param A Allocation matrix
#' @param D Dictionary matrix
#' @param model_obj Model object containing data
#'
#' @return Log-likelihood value
#' @keywords internal
#' @importFrom stats dnbinom
negative_binomial_loglikelihood_matrix <- function(A, D, model_obj) {
  start_timer("negative_binomial_loglikelihood_matrix")
  
  # Calculate proportions
  proportions <- (A %*% D) / rowSums(A)
  
  # Negative binomial likelihood
  # Y ~ NegativeBinomial(R, 1/(1 + proportions))
  loglik <- sum(stats::dnbinom(model_obj$y, size = model_obj$r, prob = 1 / (1 + proportions), log = TRUE), na.rm = TRUE)
  
  end_timer("negative_binomial_loglikelihood_matrix")
  return(loglik)
}

#' Log-likelihood for negative binomial model (vector version)
#'
#' @param propvec Proportion vector
#' @param yvec Observed counts vector
#' @param rvec Total counts vector
#'
#' @return Log-likelihood value
#' @keywords internal
#' @importFrom stats dnbinom
negative_binomial_loglikelihood_vector <- function(propvec, yvec, rvec) {
  # Negative binomial likelihood for vectors
  loglik <- sum(stats::dnbinom(yvec, size = rvec, prob = 1 / (1 + propvec), log = TRUE), na.rm = TRUE)
  return(loglik)
}

#' Initialize state for negative binomial model
#'
#' @param model_obj Model object
#' @param threshold Threshold for identifying single infections
#'
#' @return Initialized state
#' @keywords internal
#' @importFrom stats runif
negative_binomial_initialize_state <- function(model_obj, threshold = 0.001) {
  y <- model_obj$y
  r <- model_obj$r
  N <- model_obj$N
  P <- model_obj$P
  
  # Calculate ratios
  ratios <- y / r
  
  # Create initial dictionary from ratios
  cate <- (y / r > 0.5)
  if (any(is.na(cate))) {
    cate[which(is.na(cate))] <- 0
  }
  
  # Remove duplicates to create dictionary
  ratios2dict <- remove_duplicates(cate)
  assignments <- ratios2dict$assignments
  nstrain <- length(unique(assignments))
  
  # Initialize state
  state <- list()
  state$D <- ratios2dict$D
  state$A <- matrix(0, nrow = N, ncol = nstrain)
  
  # Identify single vs mixed infections
  is_single <- apply(ratios, 1, function(x) all(pmin(x, 1 - x) < threshold))
  nsingleppl <- sum(is_single, na.rm = TRUE)
  
  if (nsingleppl >= 1) {
    # Handle single infections
    which_single <- which(is_single)
    which_mixed <- c(which(!is_single), which(is.na(is_single)))
    nsinglestrain <- length(unique(assignments[which_single]))
    
    # Assign mixed infections to all strains
    state$A[which_mixed, ] <- 1
    
    # Assign single infections to specific strains
    for (isingle in 1:nsingleppl) {
      state$A[which_single[isingle], assignments[which_single[isingle]]] <- 1
    }
    
    # Reorder strains (single infection strains first)
    ord <- rep(NA, nstrain)
    ord[1:nsinglestrain] <- unique(assignments[which_single])
    ord[(nsinglestrain + 1):nstrain] <- setdiff(unique(assignments), unique(assignments[which_single]))
    
    state$A <- state$A[, ord, drop = FALSE]
    state$D <- state$D[ord, , drop = FALSE]
    state$mixed <- which_mixed
    state$kmin <- nsinglestrain + 1
    
    # Resolve exceptions
    state$loglik <- negative_binomial_loglikelihood_matrix(state$A, state$D, model_obj)
    iter <- 1
    while (is.infinite(state$loglik)) {
      state <- negative_binomial_resolve_exceptions(state, model_obj)
      iter <- iter + 1
      if (iter > 10) {
        stop("Failed to initialize valid state after 10 iterations")
      }
      state$loglik <- negative_binomial_loglikelihood_matrix(state$A, state$D, model_obj)
    }
  } else {
    # All mixed infections
    state$mixed <- 1:N
    state$kmin <- 1
    state$A[state$mixed, ] <- 1
    state$loglik <- negative_binomial_loglikelihood_matrix(state$A, state$D, model_obj)
    
    # Resolve exceptions
    iter <- 1
    while (is.infinite(state$loglik)) {
      state <- negative_binomial_resolve_exceptions(state, model_obj)
      iter <- iter + 1
      if (iter > 10) {
        stop("Failed to initialize valid state after 10 iterations")
      }
    }
  }
  
  return(state)
}

#' Resolve exceptions for negative binomial model
#'
#' @param state Current state
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
negative_binomial_resolve_exceptions <- function(state, model_obj) {
  y <- model_obj$y
  r <- model_obj$r
  
  # Calculate current proportions
  ratios <- (state$A %*% state$D) / rowSums(state$A)
  
  # When y > 0, we must have ratios > 0
  exceptions1 <- which(y != 0 & ratios == 0, arr.ind = TRUE)
  if (nrow(exceptions1) > 0) {
    # Ensure we have enough columns in the matrices
    needed_cols <- max(state$kmin + 1, ncol(state$A))
    if (needed_cols > ncol(state$A)) {
      # Expand matrices
      state$A <- cbind(state$A, matrix(0, nrow = nrow(state$A), ncol = needed_cols - ncol(state$A)))
      state$D <- rbind(state$D, matrix(0, ncol = ncol(state$D), nrow = needed_cols - nrow(state$D)))
    }
    
    state$A[exceptions1[, 1], state$kmin + 1] <- 1
    state$D[state$kmin + 1, exceptions1[, 2]] <- 1
  }
  
  # Update log-likelihood
  state$loglik <- negative_binomial_loglikelihood_matrix(state$A, state$D, model_obj)
  
  return(state)
}
