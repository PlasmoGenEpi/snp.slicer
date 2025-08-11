#' Categorical Model for SNP-Slice
#'
#' @description
#' Implementation of the categorical observation model for SNP-Slice.
#' This model handles categorical observations (0, 0.5, 1) with error parameters.

#' Log-likelihood for categorical model (matrix version)
#'
#' @param A Allocation matrix
#' @param D Dictionary matrix
#' @param model_obj Model object containing data
#'
#' @return Log-likelihood value
#' @keywords internal
categorical_loglikelihood_matrix <- function(A, D, model_obj) {
  start_timer("categorical_loglikelihood_matrix")
  
  # Calculate proportions
  proportions <- (A %*% D) / rowSums(A)
  
  # Categorical likelihood using likelihood table
  loglik <- sum(sapply(1:model_obj$P, function(p) {
    categorical_loglikelihood_vector(proportions[, p], model_obj$y[, p])
  }))
  
  end_timer("categorical_loglikelihood_matrix")
  return(loglik)
}

#' Log-likelihood for categorical model (vector version)
#'
#' @param propvec Proportion vector
#' @param yvec Observed categorical values
#' @param rvec Not used for categorical model (kept for interface consistency)
#'
#' @return Log-likelihood value
#' @keywords internal
categorical_loglikelihood_vector <- function(propvec, yvec, rvec = NULL) {
  start_timer("categorical_loglikelihood_vector")
  
  # Create likelihood table based on error parameters
  e1 <- 0.05  # Will be passed from model object in full implementation
  e2 <- 0.05
  
  logf0 <- log(c(1 - e1, e1, 0))
  logf1 <- log(c(0, e1, 1 - e1))
  logfmix <- log(c(e2/2, 1 - e2, e2/2))
  
  llik_tab <- matrix(0, nrow = 3, ncol = 3)
  llik_tab[1, ] <- logf0
  llik_tab[2, ] <- logfmix
  llik_tab[3, ] <- logf1
  
  # Convert proportions to categorical
  propvec_cat <- cut(propvec, c(-Inf, 0, 0.99, 1.1), labels = c(0, 0.5, 1))
  yvec_cat <- factor(yvec, levels = c(0, 0.5, 1))
  
  # Calculate log-likelihood
  loglik <- sum(table(propvec_cat, yvec_cat) * llik_tab, na.rm = TRUE)
  
  end_timer("categorical_loglikelihood_vector")
  return(loglik)
}

#' Initialize state for categorical model
#'
#' @param model_obj Model object
#' @param threshold Threshold for identifying single infections
#'
#' @return Initialized state
#' @keywords internal
categorical_initialize_state <- function(model_obj, threshold = 0.001) {
  y <- model_obj$y
  N <- model_obj$N
  P <- model_obj$P
  
  # Replace NA or 0.5 in y with 0 or 1 (randomly)
  cate <- y
  for (i in 1:N) {
    for (j in 1:P) {
      if (is.na(y[i, j]) || y[i, j] == 0.5) {
        cate[i, j] <- (stats::runif(1) < 0.5)
      }
    }
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
  is_single <- apply(y, 1, function(x) all(x == 0 | x == 1))
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
    remaining_strains <- setdiff(unique(assignments), unique(assignments[which_single]))
    if (length(remaining_strains) > 0) {
      ord[(nsinglestrain + 1):(nsinglestrain + length(remaining_strains))] <- remaining_strains
    }
    
    state$A <- state$A[, ord, drop = FALSE]
    state$D <- state$D[ord, , drop = FALSE]
    state$mixed <- which_mixed
    state$kmin <- nsinglestrain + 1
    
    # Resolve exceptions
    state$loglik <- categorical_loglikelihood_matrix(state$A, state$D, model_obj)
    iter <- 1
    while (is.infinite(state$loglik)) {
      state <- categorical_resolve_exceptions(state, model_obj)
      iter <- iter + 1
      if (iter > 10) {
        stop("Failed to initialize valid state after 10 iterations")
      }
      state$loglik <- categorical_loglikelihood_matrix(state$A, state$D, model_obj)
    }
  } else {
    # All mixed infections
    state$mixed <- 1:N
    state$kmin <- 1
    state$A[state$mixed, ] <- 1
    state$loglik <- categorical_loglikelihood_matrix(state$A, state$D, model_obj)
    
    # Resolve exceptions
    iter <- 1
    while (is.infinite(state$loglik)) {
      state <- categorical_resolve_exceptions(state, model_obj)
      iter <- iter + 1
      if (iter > 10) {
        stop("Failed to initialize valid state after 10 iterations")
      }
    }
  }
  
  return(state)
}

#' Resolve exceptions for categorical model
#'
#' @param state Current state
#' @param model_obj Model object
#'
#' @return Updated state
#' @keywords internal
categorical_resolve_exceptions <- function(state, model_obj) {
  y <- model_obj$y
  
  # Calculate current proportions
  ratios <- (state$A %*% state$D) / rowSums(state$A)
  
  # When y == 0, cannot have proportions == 1
  exceptions2 <- which(ratios == 1 & y == 0, arr.ind = TRUE)
  if (nrow(exceptions2) > 0) {
    state$A[exceptions2[, 1], state$kmin] <- 1
    state$D[state$kmin, exceptions2[, 2]] <- 0
  }
  
  # Recalculate ratios
  ratios <- (state$A %*% state$D) / rowSums(state$A)
  
  # When y == 1, cannot have proportions == 0
  exceptions1 <- which(y == 1 & ratios == 0, arr.ind = TRUE)
  if (nrow(exceptions1) > 0) {
    state$A[exceptions1[, 1], state$kmin + 1] <- 1
    state$D[state$kmin + 1, exceptions1[, 2]] <- 1
  }
  
  # Update log-likelihood
  state$loglik <- categorical_loglikelihood_matrix(state$A, state$D, model_obj)
  
  return(state)
}
